const labModelLabels = {
  QS: "Charge sharing (QS)",
  QR: "Charge redistribution (QR)",
  DIMC: "Digital IMC (DIMC)",
  IS: "Current summing (IS)",
  "QS-QR": "Hybrid charge sharing / redistribution (QS–QR)",
  "": "Not classified in the dataset"
};

const labModelOrder = ["QR", "QS", "DIMC", "IS", "QS-QR", ""];

function labNumber(value) {
  if (value === null || value === undefined || String(value).trim() === "") return null;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

function labEqual(left, right, tolerance = 1e-8) {
  return left !== null && right !== null && Math.abs(left - right) <= tolerance;
}

function labFormat(value, digits = 2) {
  if (!Number.isFinite(value)) return "—";
  return value.toLocaleString(undefined, { maximumFractionDigits: digits });
}

function labUniquePapers(rows) {
  return [...new Set(rows.map((row) => row.Index).filter(Boolean))];
}

function labRange(values) {
  const finite = values.filter(Number.isFinite);
  return finite.length ? [Math.min(...finite), Math.max(...finite)] : null;
}

function labClamped(target, selected) {
  return Math.abs(target - selected) > 1e-6 * Math.max(1, Math.abs(selected));
}

function labBracket(values, selected, logarithmic = false) {
  const transform = logarithmic ? Math.log : (value) => value;
  const target = Math.max(values[0], Math.min(values[values.length - 1], selected));
  let upper = values.findIndex((value) => value >= target);
  if (upper === -1) return { lower: values.length - 1, upper: values.length - 1, mix: 0, clamped: true };
  if (upper === 0) return { lower: 0, upper: 0, mix: 0, clamped: labClamped(target, selected) };
  const lower = upper - 1;
  const mix = (transform(target) - transform(values[lower])) / (transform(values[upper]) - transform(values[lower]));
  return { lower, upper, mix, clamped: labClamped(target, selected) };
}

function labBilinear(grid, rowBracket, columnBracket) {
  const r0 = grid[rowBracket.lower][columnBracket.lower];
  const r1 = grid[rowBracket.upper][columnBracket.lower];
  const c0 = grid[rowBracket.lower][columnBracket.upper];
  const c1 = grid[rowBracket.upper][columnBracket.upper];
  const lower = r0 + (r1 - r0) * rowBracket.mix;
  const upper = c0 + (c1 - c0) * rowBracket.mix;
  return lower + (upper - lower) * columnBracket.mix;
}

function labSeededRandom(seed) {
  let state = seed >>> 0;
  return () => {
    state += 0x6d2b79f5;
    let value = state;
    value = Math.imul(value ^ value >>> 15, value | 1);
    value ^= value + Math.imul(value ^ value >>> 7, value | 61);
    return ((value ^ value >>> 14) >>> 0) / 4294967296;
  };
}

function labNormal(random) {
  const first = Math.max(Number.EPSILON, random());
  const second = random();
  return Math.sqrt(-2 * Math.log(first)) * Math.cos(2 * Math.PI * second);
}

function labVariance(values) {
  const mean = values.reduce((sum, value) => sum + value, 0) / values.length;
  return values.reduce((sum, value) => sum + (value - mean) ** 2, 0) / values.length;
}

// ADC window used by the paper's sweeps: update_b_adc_mimo(mu = N / 4, B) in
// mihirvk2/tcas-mimo-imc-2025 IMC_model/imc_dp_unit.py. Limits are in units of the
// pre-ADC level spacing; the source defines them only for N = 64, 128, 256 and B = 3..9.
function labQrTunedWindow(dimension, adcBits) {
  if (![64, 128, 256].includes(dimension) || !Number.isInteger(adcBits) || adcBits < 3 || adcBits > 9) return null;
  const mu = dimension * 0.25;
  const centred = (halfWidth) => [mu - halfWidth - 0.5, mu + halfWidth - 0.5];
  if (adcBits === 3) return centred(3);
  if (adcBits === 4) return centred(7);
  return {
    5: { 64: [0.5, 30.5], 128: centred(15), 256: centred(15) },
    6: { 64: [0.5, 62.5], 128: [0.5, 62.5], 256: centred(31) },
    7: { 64: [0.25, 63.25], 128: [0.5, 126.5], 256: [0.5, 126.5] },
    8: { 64: [0.125, 63.625], 128: [0.25, 127.25], 256: [0.5, 254.5] },
    9: { 64: [0.0625, 63.8125], 128: [0.125, 127.625], 256: [0.25, 255.25] }
  }[adcBits][dimension];
}

// Uniform ADC from IMC_model/imc_adc.py: 2^B - 1 evenly spaced thresholds from t1 to tM,
// outputs at threshold + LSB/2 (below t1: t1 - LSB/2), expressed in dot-product counts.
function labQrAdc(t1, tM, adcBits, delta) {
  const count = 2 ** adcBits - 1;
  return { t1, lsb: (tM - t1) / (count - 1), count, delta };
}

function labQrConvert(adc, volts) {
  if (volts < adc.t1) return (adc.t1 - adc.lsb / 2) / adc.delta;
  const index = Math.min(adc.count - 1, Math.floor((volts - adc.t1) / adc.lsb));
  return (adc.t1 + index * adc.lsb + adc.lsb / 2) / adc.delta;
}

function labSnrDb(ideal, modeled) {
  // The source treats matches within floating-point tolerance as exact and reports 100 dB when no sample errs.
  const errors = ideal.map((value, index) => (Math.abs(value - modeled[index]) < 1e-6 ? 0 : value - modeled[index]));
  if (errors.every((error) => error === 0)) return { db: 100, errorFree: true };
  return { db: 10 * Math.log10(labVariance(ideal) / labVariance(errors)), errorFree: false };
}

function labSramQrSnr(candidate) {
  const { inputBits, weightBits, adcBits } = candidate;
  if (![inputBits, weightBits, adcBits].every(Number.isInteger) || inputBits > 16 || weightBits > 16 || adcBits < 2) return null;
  const requestedDimension = Math.max(1, Math.round(candidate.dimension));
  const dimension = Math.min(1024, requestedDimension);
  const samples = Math.max(200, Math.min(1000, Math.floor(3e7 / (dimension * inputBits * weightBits))));
  const random = labSeededRandom(280000 + dimension * 31 + inputBits * 7 + weightBits * 13 + adcBits * 17);
  const cMean = candidate.cellCap;
  const cSigma = 2.1 * 10 ** -2.5 * Math.sqrt(cMean);
  const vdd = 0.9;
  const cPar = cMean * dimension * 0.3 + 2.04278;
  const delta = cMean * vdd / (dimension * cMean + cPar);
  const levels = 2 ** adcBits;
  const fullAdc = labQrAdc(0.5 * delta * dimension / levels, (levels - 1.5) * delta * dimension / levels, adcBits, delta);
  const window = labQrTunedWindow(dimension, adcBits);
  const tunedAdc = window ? labQrAdc(window[0] * delta, window[1] * delta, adcBits, delta) : null;
  const capacitance = new Float64Array(dimension * weightBits);
  for (let index = 0; index < capacitance.length; index += 1) capacitance[index] = Math.max(1e-6, cMean + cSigma * labNormal(random));
  const capTotals = Array.from({ length: weightBits }, (_, bit) => {
    let total = cPar;
    for (let row = 0; row < dimension; row += 1) total += capacitance[row * weightBits + bit];
    return total;
  });
  const noiseVolts = candidate.adcNoiseMv * 1e-3;
  const weightCodes = new Int32Array(dimension);
  const inputCodes = new Int32Array(dimension);
  const ideal = new Float64Array(samples);
  const full = new Float64Array(samples);
  const tuned = new Float64Array(samples);

  for (let sample = 0; sample < samples; sample += 1) {
    let product = 0;
    for (let row = 0; row < dimension; row += 1) {
      const weight = Math.floor(random() * 2 ** weightBits) - 2 ** (weightBits - 1);
      const input = Math.floor(random() * 2 ** inputBits) - 2 ** (inputBits - 1);
      product += weight * input;
      weightCodes[row] = weight < 0 ? weight + 2 ** weightBits : weight;
      inputCodes[row] = input < 0 ? input + 2 ** inputBits : input;
    }
    ideal[sample] = product;
    let fullTotal = 0;
    let tunedTotal = 0;
    for (let inputBit = 0; inputBit < inputBits; inputBit += 1) {
      const inputShift = inputBits - 1 - inputBit;
      let fullColumns = 0;
      let tunedColumns = 0;
      for (let weightBit = 0; weightBit < weightBits; weightBit += 1) {
        const weightShift = weightBits - 1 - weightBit;
        let activeCap = 0;
        for (let row = 0; row < dimension; row += 1) {
          if ((weightCodes[row] >> weightShift) & (inputCodes[row] >> inputShift) & 1) activeCap += capacitance[row * weightBits + weightBit];
        }
        const volts = activeCap * vdd / capTotals[weightBit] + noiseVolts * labNormal(random);
        const place = weightBit === 0 ? -(2 ** weightShift) : 2 ** weightShift;
        fullColumns += place * labQrConvert(fullAdc, volts);
        if (tunedAdc) tunedColumns += place * labQrConvert(tunedAdc, volts);
      }
      const inputPlace = inputBit === 0 ? -(2 ** inputShift) : 2 ** inputShift;
      fullTotal += inputPlace * fullColumns;
      tunedTotal += inputPlace * tunedColumns;
    }
    full[sample] = fullTotal;
    tuned[sample] = tunedTotal;
  }
  return {
    dimension, requestedDimension, samples, window,
    fullRange: labSnrDb(ideal, full),
    tunedWindow: tunedAdc ? labSnrDb(ideal, tuned) : null
  };
}

// Energy constants from Cqr_sweep.py / Badc_sweep.py: E_ADC = 0.1·B + 1e-6·4^B pJ and
// E_IA = 1.244·C_QR·N·VDD² fJ per column conversion. Normalizing one conversion over its
// 2N one-bit operations follows the benchmark's N_1b convention; that split is the lab's.
function labSramQrEnergy(candidate) {
  const adcBits = candidate.adcBits;
  if (!Number.isFinite(adcBits) || adcBits <= 0) return null;
  const vdd = 0.9;
  const dimension = Math.max(1, candidate.dimension);
  const adcFj = (0.1 * adcBits + 1e-6 * 4 ** adcBits) * 1000;
  const arrayFj = 1.244 * candidate.cellCap * dimension * vdd ** 2;
  return { adcFj, arrayFj, perOpFj: (adcFj + arrayFj) / (2 * dimension) };
}

const labAsimWorkloads = {
  "resnet18-cifar10": { label: "ResNet-18 · CIFAR-10", task: "CIFAR-10", model: "ResNet-18", chance: 10, checkpoint: "resnet18_cifar10_w8a8_pact_trim_100.pkl" },
  "resnet18-imagenet": { label: "ResNet-18 · ImageNet", task: "ImageNet", model: "ResNet-18", chance: 0.1, checkpoint: "resnet18_imagenet_w8a8_pact_trim_50.pkl" },
  "vitb32-cifar10": { label: "ViT-B-32 · CIFAR-10", task: "CIFAR-10", model: "ViT-B-32", chance: 10, checkpoint: "vitb32_cifar10_w8a8_trim_60.pkl" },
  "vitb32-imagenet": { label: "ViT-B-32 · ImageNet", task: "ImageNet", model: "ViT-B-32", chance: 0.1, checkpoint: "vitb32_imagenet_w8a8_trim_40.pkl" }
};

const labAsimEncodings = {
  0: { key: "bitSerial", label: "bit-serial" },
  2: { key: "2", label: "2-bit parallel" },
  4: { key: "4", label: "4-bit parallel" }
};

function labInterpolate(points, x) {
  const first = points[0];
  const last = points[points.length - 1];
  const onGrid = points.some(([px]) => Math.abs(px - x) < 1e-9);
  if (x <= first[0]) return { value: first[1], onGrid, beyond: x < first[0] };
  if (x >= last[0]) return { value: last[1], onGrid, beyond: x > last[0] };
  const upper = points.findIndex(([px]) => px >= x);
  const [x0, y0] = points[upper - 1];
  const [x1, y1] = points[upper];
  return { value: y0 + (y1 - y0) * (x - x0) / (x1 - x0), onGrid, beyond: false };
}

// Settings that trip assertions in the ASiM modules (asim_conv.py, asim_multiheadattention.py).
function labAsimProblems(candidate, workload) {
  const problems = [];
  if (candidate.weightBits <= 1) problems.push("1-bit weights divide by zero in asim_conv.py.");
  if (candidate.asimEncoding > candidate.inputBits) problems.push("activation encoding cannot exceed the activation precision.");
  if (workload.model === "ResNet-18" && Math.round(candidate.dimension) < 9) problems.push("convolution layers need nrow of at least 9.");
  if (workload.model === "ViT-B-32" && (candidate.inputBits <= 1 || (candidate.asimEncoding > 0 && candidate.asimEncoding >= candidate.inputBits))) problems.push("attention needs more than 1 activation bit and key encoding below the query precision.");
  return problems;
}

function labAsimConfig(candidate, workload, problems) {
  const rows = Math.round(candidate.dimension);
  const encoding = candidate.asimEncoding === 0 ? "None" : String(candidate.asimEncoding);
  const { inputBits: xbit, weightBits: wbit, adcBits: adc } = candidate;
  const sigmas = (prefix) => [
    `cfg.${prefix}_rand_noise_sigma = ${candidate.asimRandomNoise}`,
    `cfg.${prefix}_non_linear_sigma = ${candidate.asimNonlinearity}`
  ];
  const lines = [
    "# Screening translation for ASiM examples/simulation/<workload>/config.py",
    ...problems.map((problem) => `# Will not run as set: ${problem}`),
    `cfg.task = '${workload.task}'`,
    candidate.asimTraining === "nat"
      ? `cfg.pretrain_model_path = r'<ASiM_Models>/${workload.checkpoint}'  # released post-NAT W8A8 checkpoint`
      : "cfg.pretrain_model_path = r'<your QAT checkpoint>.pkl'  # ASiM releases post-NAT checkpoints only",
    "# The repository's example configs use nrow = 255."
  ];
  if (workload.model === "ResNet-18") {
    ["conv", "linear"].forEach((module) => {
      lines.push(
        `cfg.asim_${module}_wbit = ${wbit}`,
        `cfg.asim_${module}_xbit = ${xbit}`,
        `cfg.asim_${module}_adc_prec = ${adc}`,
        `cfg.asim_${module}_nrow = ${rows}`,
        ...sigmas(`asim_${module}`),
        `cfg.asim_${module === "conv" ? "act_enc" : "linear_act_enc"} = ${encoding}`,
        `cfg.asim_${module}_trim_noise = 0.0`
      );
    });
  } else {
    lines.push(
      `cfg.asim_vit_attn_nrow = ${rows}`,
      `cfg.asim_vit_attn_qk_qbit = ${xbit}`,
      `cfg.asim_vit_attn_qk_kbit = ${xbit + 1}  # sign bit plus magnitude; the example uses 9 with 8-bit activations`,
      `cfg.asim_vit_attn_av_abit = ${xbit}`,
      `cfg.asim_vit_attn_av_vbit = ${xbit}`,
      `cfg.asim_vit_attn_proj_wbit = ${wbit}`,
      `cfg.asim_vit_attn_proj_xbit = ${xbit + 1}  # the example uses 9 with 8-bit activations`,
      ...["qk", "av", "proj"].map((part) => `cfg.asim_vit_attn_${part}_adc_prec = ${adc}`),
      ...["qk", "av", "proj"].flatMap((part) => sigmas(`asim_vit_attn_${part}`)),
      `cfg.asim_vit_attn_qk_k_enc = ${encoding}`,
      `cfg.asim_vit_attn_av_a_enc = ${encoding}`,
      `cfg.asim_vit_attn_proj_act_enc = ${encoding}`,
      "cfg.asim_vit_attn_attn_trim_noise = 0.0",
      "cfg.asim_vit_attn_proj_trim_noise = 0.0"
    );
    ["mlp", "fc"].forEach((module) => {
      lines.push(
        `cfg.asim_vit_${module}_wbit = ${wbit}`,
        `cfg.asim_vit_${module}_xbit = ${xbit}`,
        `cfg.asim_vit_${module}_adc_prec = ${adc}`,
        `cfg.asim_vit_${module}_nrow = ${rows}`,
        ...sigmas(`asim_vit_${module}`),
        `cfg.asim_vit_${module}_act_enc = ${encoding}`,
        `cfg.asim_vit_${module}_trim_noise = 0.0`
      );
    });
    lines.push(`cfg.asim_vit_quant_conv_wbit = ${wbit}`, `cfg.asim_vit_quant_conv_xbit = ${xbit}`);
  }
  lines.push("# Run main/src_simulation.py to obtain the validation-set Test Acc.");
  return lines.join("\n");
}

// Reads ASiM's published tables (window.IMC_ASIM_DATA, from scripts/build-asim-data.py).
// ADC precision comes from Figs. 6-7 (no analog noise). Noise enters as the fraction of the
// above-chance margin kept in Fig. 9a/9b (or Fig. 18 without NAT), which were all measured at
// an 8-bit ADC. Anything the tables do not cover is translated or reported as an upper bound.
function labAsimScreen(candidate) {
  const workload = labAsimWorkloads[candidate.asimWorkload];
  const data = window.IMC_ASIM_DATA?.workloads?.[candidate.asimWorkload];
  if (!workload || !data) return null;
  const encoding = labAsimEncodings[candidate.asimEncoding];
  const chance = workload.chance;
  const notes = [];
  let bound = false;
  let translated = false;

  const rowBits = Math.ceil(Math.log2(Math.max(2, candidate.dimension)));
  const boundaryAdc = rowBits + candidate.asimEncoding;
  const equivalentAdc = candidate.adcBits + 8 - rowBits;
  const adcKey = Math.floor(equivalentAdc + 1e-9);
  if (rowBits !== 8) {
    translated = true;
    notes.push(`N = ${Math.round(candidate.dimension)} is mapped onto the paper's 256-row sweep by ADC headroom above the ⌈log₂N⌉-bit boundary, so ${labFormat(candidate.adcBits, 2)} bits reads as ${labFormat(equivalentAdc, 2)} bits.`);
  }
  if (Math.abs(adcKey - equivalentAdc) > 1e-9) {
    translated = true;
    notes.push("Fractional ADC precision is rounded down.");
  }
  const adcPoints = data.adc[encoding.key];
  const adcLow = adcPoints[0][0];
  const adcHigh = adcPoints[adcPoints.length - 1][0];
  let adcTop1;
  if (adcKey < adcLow) {
    adcTop1 = adcPoints[0][1];
    bound = true;
    notes.push(`That is below the tabulated ${adcLow}–${adcHigh}-bit range for ${encoding.label} encoding, so accuracy can be no better than the ${adcLow}-bit result.`);
  } else if (adcKey > adcHigh) {
    adcTop1 = adcPoints[adcPoints.length - 1][1];
    notes.push(`That is above the tabulated range, so the ${adcHigh}-bit result is used.`);
  } else {
    adcTop1 = adcPoints.find(([bits]) => bits === adcKey)[1];
  }

  const nat = candidate.asimTraining === "nat";
  const randomTable = nat || !data.nonNatRandomNoise
    ? { points: data.randomNoise[encoding.key], figure: "Fig. 9a", natData: true }
    : { points: data.nonNatRandomNoise[encoding.key], figure: "Fig. 18", natData: false };
  const noise = [
    { label: "Random noise", sigma: candidate.asimRandomNoise, ...randomTable },
    { label: "Nonlinearity", sigma: candidate.asimNonlinearity, points: data.nonlinearity[encoding.key], figure: "Fig. 9b", natData: true }
  ].map((source) => {
    if (source.sigma <= 0) return { ...source, retention: 1, top1: null };
    const reading = labInterpolate(source.points, source.sigma);
    const clean = source.points[0][1];
    if (!reading.onGrid) translated = true;
    if (reading.beyond) {
      bound = true;
      notes.push(`${source.label} above 0.2% Vpp is outside ${source.figure}; its 0.2% result is an upper bound.`);
    }
    if (!nat && source.natData) {
      bound = true;
      notes.push(`${source.label} has no table without NAT for this workload, and NAT results overstate the robustness of a QAT-only checkpoint.`);
    }
    if (clean - chance < 1) {
      bound = true;
      notes.push(`${source.label}: at an 8-bit ADC the ${encoding.label} column of ${source.figure} is already at chance, so its noise effect cannot be separated from ADC loss.`);
      return { ...source, retention: 1, top1: reading.value };
    }
    return { ...source, retention: Math.max(0, Math.min(1, (reading.value - chance) / (clean - chance))), top1: reading.value };
  });
  const noisy = noise.filter((source) => source.sigma > 0);
  if (noisy.length && equivalentAdc > 8) {
    translated = true;
    notes.push("Noise tables were measured at an 8-bit ADC; the paper reports that higher precision reduces noise damage, so this noise loss is conservative.");
  } else if (noisy.length && equivalentAdc < 8) {
    bound = true;
    notes.push("Noise tables were measured at an 8-bit ADC; lower precision increases noise damage, so the result is an upper bound.");
  }
  if (noisy.length === 2) {
    bound = true;
    notes.push("Random noise and nonlinearity were swept separately; together they can do no better than the more damaging one.");
  }
  const outside = candidate.weightBits !== 8 || candidate.inputBits !== 8;
  if (outside) notes.push("Every table uses W8A8 checkpoints; other precisions need a retrained checkpoint and a direct ASiM run.");

  const scale = (top1) => chance + (candidate.asimBaseline - chance) * Math.max(0, top1 - chance) / (data.baseline - chance);
  const adcOnly = Math.max(chance, Math.min(candidate.asimBaseline, scale(adcTop1)));
  const retention = Math.min(...noise.map((source) => source.retention));
  return {
    workload, encoding, boundaryAdc, equivalentAdc, adcOnly, bound, notes,
    paperBaseline: data.baseline,
    estimate: chance + (adcOnly - chance) * retention,
    noise: noise.map((source) => ({ ...source, scaled: source.top1 === null ? null : Math.min(candidate.asimBaseline, scale(source.top1)) })),
    evidenceClass: outside ? "Outside the tabulated W8A8 setting" : bound ? "Upper bound from ASiM tables" : translated ? "Translated from ASiM tables" : "Tabulated ASiM point",
    config: labAsimConfig(candidate, workload, labAsimProblems(candidate, workload))
  };
}

// SNDR at a non-6-bit ADC: the extra quantization noise seen in the source's N = 64 ADC sweep
// (relative to its own 6-bit point) is added as noise power to the 6-bit surface value.
function labEnvmAdcSndr(source, device, vbl, sndrSixBit, adcBits) {
  if (Math.abs(adcBits - 6) < 1e-9) return { db: sndrSixBit, method: "surface" };
  const sweep = source.adcSweep;
  const grid = sweep?.sndrDb?.[device];
  if (!grid) return { db: null, reason: `the source has no usable ADC sweep for ${device}` };
  const bits = sweep.adcBits;
  if (adcBits < bits[0] || adcBits > bits[bits.length - 1]) return { db: null, reason: `the ADC sweep covers ${bits[0]}–${bits[bits.length - 1]} bits` };
  const vblBracket = labBracket(sweep.vbl, vbl, true);
  const selected = labBilinear(grid, labBracket(bits, adcBits), vblBracket);
  const sixBit = labBilinear(grid, labBracket(bits, 6), vblBracket);
  const extraNoise = Math.max(0, 10 ** (-selected / 10) - 10 ** (-sixBit / 10));
  return { db: -10 * Math.log10(10 ** (-sndrSixBit / 10) + extraNoise), method: "adc-sweep" };
}

// Current-summing eNVM column from calmyor/eNVM-IMC-Modeling SNDRd_vs_Energy_plot.py:
// VDD·VBL·T·Geq + (k1·B + k2·4^B) per conversion, normalized over 2N one-bit operations.
function labEnvmModel(candidate) {
  const source = window.IMC_MODEL_DATA;
  const device = source?.devices?.[candidate.envmDevice];
  if (!device) return null;
  const dimensions = source.dimensions;
  const rowBracket = labBracket(dimensions, candidate.dimension, true);
  const columnBracket = labBracket(device.vbl, candidate.vbl, true);
  const modelDimension = Math.max(dimensions[0], Math.min(dimensions[dimensions.length - 1], candidate.dimension));
  const geq = labBilinear(device.averageGeq, rowBracket, columnBracket);
  const sndrSixBit = labBilinear(device.sndrDb, rowBracket, columnBracket);
  const adcBits = Number.isFinite(candidate.adcBits) && candidate.adcBits > 0 ? candidate.adcBits : 6;
  const constants = source.energyConstants;
  const periodSeconds = constants.sourceCorePeriodSeconds;
  const arrayFj = device.vdd * candidate.vbl * periodSeconds * geq * 1e15;
  const adcFj = (constants.adcLinearJPerBit * adcBits + constants.adcExponentialJ * 4 ** adcBits) * 1e15;
  return {
    modelDimension, dimensionClamped: rowBracket.clamped, vblClamped: columnBracket.clamped,
    adcBits, assumedAdc: adcBits !== candidate.adcBits, periodNs: periodSeconds * 1e9,
    arrayFj, adcFj, perOpFj: (arrayFj + adcFj) / (2 * modelDimension),
    sndrSixBit, sndr: labEnvmAdcSndr(source, candidate.envmDevice, candidate.vbl, sndrSixBit, adcBits)
  };
}

function labModelEnergyPerOp(candidate) {
  if (candidate.architecture === "eNVM" && candidate.model === "IS") return labEnvmModel(candidate)?.perOpFj ?? null;
  if (candidate.architecture === "SRAM" && candidate.model === "QR") return labSramQrEnergy(candidate)?.perOpFj ?? null;
  return null;
}

let labQrCache = { key: "", result: null };
function labSramQrSnrCached(candidate) {
  const key = [candidate.dimension, candidate.inputBits, candidate.weightBits, candidate.adcBits, candidate.cellCap, candidate.adcNoiseMv].join("|");
  if (labQrCache.key !== key) labQrCache = { key, result: labSramQrSnr(candidate) };
  return labQrCache.result;
}

async function initDetailedLab() {
  const chart = document.querySelector("#lab-chart");
  if (!chart) return;
  setupReportedPointTooltip(chart);
  const rawRows = await loadBenchmarkRows();
  const rows = rawRows.map((row) => ({
    ...row,
    Architecture: row.Architecture.trim(),
    model: row["Compute Model"].trim(),
    tech: labNumber(row["Tech (nm)"]),
    dimension: labNumber(row.N),
    inputBits: labNumber(row.B_x),
    weightBits: labNumber(row.B_w),
    adcBits: labNumber(row.B_ADC) || null,
    efficiency: labNumber(row["TOPS/W"]),
    density: labNumber(row["TOPS/mm2"]),
    throughput: labNumber(row.TOPS)
  }));
  const plottedRows = rows.filter((row) => row.efficiency > 0 && row.density > 0);

  const controls = {
    architecture: document.querySelector("#lab-architecture"), model: document.querySelector("#lab-model"),
    tech: document.querySelector("#lab-tech"), dimension: document.querySelector("#lab-dimension"),
    inputBits: document.querySelector("#lab-input-bits"), weightBits: document.querySelector("#lab-weight-bits"),
    adcBits: document.querySelector("#lab-adc-bits"), adcColumns: document.querySelector("#lab-adc-columns"),
    coreLatency: document.querySelector("#lab-core-latency"), area: document.querySelector("#lab-area"),
    energySource: document.querySelector("#lab-energy-source"), columnEnergy: document.querySelector("#lab-column-energy"),
    information: document.querySelector("#lab-information"),
    hasSnr: document.querySelector("#lab-has-snr"), snr: document.querySelector("#lab-snr"),
    envmDevice: document.querySelector("#lab-envm-device"), vbl: document.querySelector("#lab-vbl"),
    cellCap: document.querySelector("#lab-cell-cap"), adcNoise: document.querySelector("#lab-adc-noise"),
    asimWorkload: document.querySelector("#lab-asim-workload"), asimBaseline: document.querySelector("#lab-asim-baseline"),
    asimEncoding: document.querySelector("#lab-asim-encoding"), asimRandomNoise: document.querySelector("#lab-asim-random-noise"),
    asimNonlinearity: document.querySelector("#lab-asim-nonlinearity"), asimTraining: document.querySelector("#lab-asim-training")
  };
  const outputs = {
    tech: document.querySelector("#lab-tech-output"), coreLatency: document.querySelector("#lab-core-latency-output"),
    area: document.querySelector("#lab-area-output"), columnEnergy: document.querySelector("#lab-column-energy-output"),
    information: document.querySelector("#lab-information-output"), snr: document.querySelector("#lab-snr-output"),
    vbl: document.querySelector("#lab-vbl-output"), cellCap: document.querySelector("#lab-cell-cap-output"),
    adcNoise: document.querySelector("#lab-adc-noise-output"),
    asimRandomNoise: document.querySelector("#lab-asim-random-noise-output"), asimNonlinearity: document.querySelector("#lab-asim-nonlinearity-output")
  };
  const colors = { SRAM: "#218a5b", eNVM: "#df4b45", eFlash: "#df4b45", Digital: "#3268cc", eDRAM: "#d98a1c" };
  const width = Math.max(320, Math.round(chart.getBoundingClientRect().width || 900));
  const compact = width < 560;
  const height = compact ? 420 : 560;
  const margin = compact ? { top: 20, right: 16, bottom: 66, left: 72 } : { top: 24, right: 28, bottom: 78, left: 96 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const xDomain = [0.003, 6000];
  const yDomain = [4, 25000];
  const xTicks = compact ? [0.01, 1, 100, 1000] : [0.01, 0.1, 1, 10, 100, 1000];
  const yTicks = [10, 100, 1000, 10000];
  const position = (value, domain, size) => (Math.log10(value) - Math.log10(domain[0])) / (Math.log10(domain[1]) - Math.log10(domain[0])) * size;
  const xPosition = (value) => margin.left + position(value, xDomain, plotWidth);
  const yPosition = (value) => margin.top + plotHeight - position(value, yDomain, plotHeight);
  chart.setAttribute("viewBox", `0 0 ${width} ${height}`);

  function positiveNumber(control, fallback) {
    const value = Number(control.value);
    return control.value.trim() !== "" && Number.isFinite(value) && value > 0 ? value : fallback;
  }

  // Placement follows the benchmark CSV's own normalization:
  // N_1b = 2·N_ADC·N·B_x·B_w, TOPS = N_1b / T_core, E_OP1 = E_col / (2·N·B_x·B_w).
  function candidateValues() {
    const candidate = {
      architecture: controls.architecture.value, model: controls.model.value,
      tech: Number(controls.tech.value), dimension: Math.round(positiveNumber(controls.dimension, 1)),
      inputBits: positiveNumber(controls.inputBits, 1), weightBits: positiveNumber(controls.weightBits, 1),
      adcBits: controls.adcBits.disabled || controls.adcBits.value.trim() === "" ? null : positiveNumber(controls.adcBits, null),
      adcColumns: Math.round(positiveNumber(controls.adcColumns, 1)),
      coreLatencyNs: 10 ** Number(controls.coreLatency.value), areaMm2: 10 ** Number(controls.area.value),
      enteredColumnEnergyFj: 10 ** Number(controls.columnEnergy.value), energySource: controls.energySource.value,
      information: Number(controls.information.value), hasSnr: controls.hasSnr.checked, snr: Number(controls.snr.value),
      envmDevice: controls.envmDevice.value, vbl: 10 ** Number(controls.vbl.value),
      cellCap: 10 ** Number(controls.cellCap.value), adcNoiseMv: Number(controls.adcNoise.value),
      asimWorkload: controls.asimWorkload.value, asimBaseline: Number(controls.asimBaseline.value),
      asimEncoding: Number(controls.asimEncoding.value), asimRandomNoise: Number(controls.asimRandomNoise.value),
      asimNonlinearity: Number(controls.asimNonlinearity.value), asimTraining: controls.asimTraining.value
    };
    const opsPerColumn = 2 * candidate.dimension * candidate.inputBits * candidate.weightBits;
    candidate.bitOps = candidate.adcColumns * opsPerColumn;
    candidate.modelEnergyPerOpFj = labModelEnergyPerOp(candidate);
    candidate.usesModelEnergy = candidate.energySource === "model" && Number.isFinite(candidate.modelEnergyPerOpFj);
    candidate.energyPerOpFj = candidate.usesModelEnergy ? candidate.modelEnergyPerOpFj : candidate.enteredColumnEnergyFj / opsPerColumn;
    candidate.columnEnergyFj = candidate.energyPerOpFj * opsPerColumn;
    candidate.throughput = candidate.bitOps / (candidate.coreLatencyNs * 1e-9) / 1e12;
    candidate.efficiency = 1000 / candidate.energyPerOpFj;
    candidate.density = candidate.throughput / candidate.areaMm2;
    return candidate;
  }

  function refreshModelOptions() {
    const architecture = controls.architecture.value;
    const counts = new Map();
    rows.filter((row) => row.Architecture === architecture).forEach((row) => counts.set(row.model, (counts.get(row.model) || 0) + 1));
    if (architecture === "eNVM" && !counts.has("IS")) counts.set("IS", 0);
    const prior = controls.model.options.length ? controls.model.value : null;
    controls.model.replaceChildren(...labModelOrder.filter((model) => counts.has(model)).map((model) => new Option(labModelLabels[model] || model, model)));
    // Keep a classified model the new architecture also reports; otherwise open on its most-reported model.
    controls.model.value = prior && counts.has(prior) ? prior : [...counts.entries()].sort((a, b) => b[1] - a[1])[0][0];
    refreshModelScope();
  }

  function refreshModelScope() {
    const architecture = controls.architecture.value;
    const model = controls.model.value;
    const modelRows = rows.filter((row) => row.Architecture === architecture && row.model === model);
    const dimensions = [...new Set(modelRows.map((row) => row.dimension).filter(Number.isFinite))].sort((a, b) => a - b);
    document.querySelector("#lab-dimension-help").textContent = dimensions.length ? `Reported for this pair: ${dimensions.join(", ")}.` : "No reported dimension for this pair.";
    const noAdc = model === "DIMC" || architecture === "Digital";
    const wasDisabled = controls.adcBits.disabled;
    controls.adcBits.disabled = noAdc;
    if (noAdc) controls.adcBits.value = "";
    else if (wasDisabled && !controls.adcBits.value) controls.adcBits.value = modelRows.some((row) => row.adcBits === 6) ? "6" : "5";
    document.querySelector("#envm-model-controls").hidden = !(architecture === "eNVM" && model === "IS");
    document.querySelector("#sram-model-controls").hidden = !(architecture === "SRAM" && model === "QR");
    document.querySelector("#asim-model-controls").hidden = !(architecture === "SRAM" && ["QS", "QR", "QS-QR"].includes(model));
  }

  function updateEvidence(candidate) {
    const pair = rows.filter((row) => row.Architecture === candidate.architecture && row.model === candidate.model);
    const node = pair.filter((row) => labEqual(row.tech, candidate.tech));
    const dimension = node.filter((row) => labEqual(row.dimension, candidate.dimension));
    const precision = dimension.filter((row) => labEqual(row.inputBits, candidate.inputBits) && labEqual(row.weightBits, candidate.weightBits));
    const sameAdc = (row) => (candidate.adcBits === null ? row.adcBits === null : labEqual(row.adcBits, candidate.adcBits));
    const exact = precision.filter(sameAdc);
    const sameExceptNode = pair.filter((row) => labEqual(row.dimension, candidate.dimension) && labEqual(row.inputBits, candidate.inputBits) && labEqual(row.weightBits, candidate.weightBits) && sameAdc(row));
    const status = document.querySelector("#evidence-status");
    const trail = `${pair.length} model point${pair.length === 1 ? "" : "s"} → ${node.length} at ${candidate.tech} nm → ${dimension.length} at N=${candidate.dimension} → ${precision.length} at ${labFormat(candidate.inputBits, 2)}b×${labFormat(candidate.weightBits, 2)}b → ${exact.length} exact`;
    if (exact.length) {
      const papers = labUniquePapers(exact).map((index) => `#${index}`).join(", ");
      status.className = "evidence-status supported";
      status.innerHTML = `<strong>Reported configuration family</strong><p>${escapeHtml(trail)}. Matching paper${labUniquePapers(exact).length === 1 ? "" : "s"}: ${escapeHtml(papers)}.</p>`;
    } else if (sameExceptNode.length) {
      const nodes = [...new Set(sameExceptNode.map((row) => row.tech))].sort((a, b) => a - b).join(", ");
      const papers = labUniquePapers(sameExceptNode).map((index) => `#${index}`).join(", ");
      status.className = "evidence-status caution";
      status.innerHTML = `<strong>New technology-node extrapolation</strong><p>All other selected fields match ${escapeHtml(papers)} at ${escapeHtml(nodes)} nm, but not at ${candidate.tech} nm. Device and interconnect nonidealities may increase, so compute SNR can drop; the lab does not scale SNR automatically.</p><small>${escapeHtml(trail)}</small>`;
    } else {
      const closestStage = precision.length ? "ADC configuration" : dimension.length ? "precision" : node.length ? "dot-product dimension" : "technology node";
      status.className = "evidence-status caution";
      status.innerHTML = `<strong>Unreported combination — treat as a new design</strong><p>The benchmark has no paper matching the complete selection. The first unsupported step in this path is ${escapeHtml(closestStage)}; performance and accuracy require a new circuit/model evaluation.</p><small>${escapeHtml(trail)}</small>`;
    }
    return { pair, exact, sameExceptNode };
  }

  function energyComparison(candidate, modelPerOpFj, name) {
    if (candidate.usesModelEnergy) return `Placement uses the ${name} energy estimate.`;
    const ratio = candidate.energyPerOpFj / modelPerOpFj;
    return `Entered energy is ${labFormat(ratio >= 1 ? ratio : 1 / ratio, 2)}× ${ratio >= 1 ? "above" : "below"} the ${name} estimate.`;
  }

  function validateModel(candidate) {
    if (candidate.architecture === "eNVM" && candidate.model === "IS" && window.IMC_MODEL_DATA) {
      const model = labEnvmModel(candidate);
      const scope = [
        model.dimensionClamped ? `N = ${labFormat(model.modelDimension, 0)} (the edge of the model's 16–5,000 range; ${candidate.dimension} selected)` : `N = ${candidate.dimension}`,
        `${labFormat(candidate.vbl, 4)} V bitline${model.vblClamped ? " (clamped to the published range)" : ""}`,
        `the source's fixed ${labFormat(model.periodNs, 0)} ns read period`,
        `a ${labFormat(model.adcBits, 2)}-bit ADC${model.assumedAdc ? " (assumed; none selected)" : ""}`
      ].join(", ");
      const sndrText = model.sndr.db === null ? `Not modeled: ${model.sndr.reason}` : `${labFormat(model.sndr.db, 1)} dB${model.sndr.method === "adc-sweep" ? " (6-bit surface adjusted with the N = 64 ADC sweep)" : ""}`;
      return {
        lead: energyComparison(candidate, model.perOpFj, "eNVM model"),
        detail: `${candidate.envmDevice} current-summing column at ${scope}. Energy covers the column current and the ADC only; drivers, bias, and other peripherals are excluded. Behavioral model calibrated against a 22 nm MRAM prototype, with no technology-node scaling.`,
        metrics: [
          ["Model energy per 1b op", `${labFormat(model.perOpFj, 3)} fJ`],
          ["Array / ADC energy", `${labFormat(model.arrayFj, 1)} / ${labFormat(model.adcFj, 1)} fJ per conversion`],
          ["SNDR, 6-bit surface", `${labFormat(model.sndrSixBit, 1)} dB`],
          [`SNDR at ${labFormat(model.adcBits, 2)} bits`, sndrText]
        ],
        modeledSnr: model.sndr.db,
        referenceSnr: model.sndrSixBit,
        referenceNode: "22 nm MRAM calibration"
      };
    }
    if (candidate.architecture === "SRAM" && candidate.model === "QR") {
      const snr = labSramQrSnrCached(candidate);
      const energy = labSramQrEnergy(candidate);
      const atSourceNode = candidate.tech === 28;
      const snrText = (result) => (result.errorFree ? `No errors in ${snr.samples} samples (100 dB)` : `${labFormat(result.db, 1)} dB`);
      const metrics = [];
      if (snr) {
        metrics.push(["Compute SNR, full-range ADC", snrText(snr.fullRange)]);
        metrics.push(["Compute SNR, tuned ADC window", snr.tunedWindow ? snrText(snr.tunedWindow) : "Defined only for N = 64, 128, 256 with 3–9-bit ADCs"]);
        metrics.push(["Monte Carlo", `${snr.samples} samples${snr.dimension < snr.requestedDimension ? `, run at N = 1,024 (capped from ${snr.requestedDimension.toLocaleString()})` : ""}`]);
      }
      if (energy) metrics.push(["Model energy per 1b op", `${labFormat(energy.perOpFj, 3)} fJ (ADC ${labFormat(energy.adcFj, 0)} + array ${labFormat(energy.arrayFj, 1)} fJ per conversion)`]);
      return {
        lead: snr ? `${labFormat(snr.fullRange.db, 1)} dB full-range${snr.tunedWindow ? ` · ${labFormat(snr.tunedWindow.db, 1)} dB tuned` : ""} compute SNR` : "QR SNR needs integer input, weight, and ADC precision",
        detail: [
          snr ? "Capacitor mismatch and ADC thermal noise at the repository's 28 nm, 0.9 V parameterization, with no node scaling." : "The QR simulator decomposes integer input, weight, and ADC precisions; the selection has a fractional value or no ADC.",
          snr?.tunedWindow ? "The tuned window centres the ADC on the mean count N/4, as in the paper's MIMO sweeps; rare clipping errors make it vary by a few dB between runs." : "",
          energy ? `${energyComparison(candidate, energy.perOpFj, "QR model")} Splitting the repository's ADC and array constants over 2N one-bit operations is the lab's normalization.` : "",
          snr && !atSourceNode ? `The SNR is not transferred to ${candidate.tech} nm.` : ""
        ].filter(Boolean).join(" "),
        metrics,
        modeledSnr: snr && atSourceNode ? snr.fullRange.db : null,
        referenceSnr: snr ? snr.fullRange.db : null,
        referenceNode: "28 nm"
      };
    }
    if (candidate.architecture === "SRAM" && candidate.model === "DIMC") {
      const referenceLatencyUs = 2150;
      const referenceN = 1024;
      const referenceTops = 2 * referenceN ** 3 / (referenceLatencyUs * 1e-6) / 1e12;
      return {
        lead: "Behavioral DIMC reference, not a performance certificate",
        detail: "The Explorer BackdooringDiT path uses 8-bit weight/activation fake quantization. Its 2,150 µs DIMC latency is a forced simulation value extrapolated from an approximately 1-TOPS design, not a measured energy or area model.",
        metrics: [
          ["Reference operation", "1024×1024 INT8 matmul"],
          ["Forced latency", "2,150 µs"],
          ["Implied throughput", `${labFormat(referenceTops, 2)} TOPS`],
          ["Energy / density validation", "Not provided"]
        ],
        modeledSnr: null,
        referenceSnr: null
      };
    }
    return {
      lead: "Benchmark evidence only",
      detail: "No attached external model covers this architecture and compute-model pair. The lab can check reported coordinates and exact paper matches, but it cannot certify the selected energy, area, throughput, or SNR.",
      metrics: [], modeledSnr: null, referenceSnr: null
    };
  }

  function updateAccuracy(candidate, modelResult) {
    const title = document.querySelector("#accuracy-title");
    const lead = document.querySelector("#accuracy-lead");
    const detail = document.querySelector("#accuracy-detail");
    const marker = document.querySelector("#accuracy-marker");
    const scale = document.querySelector("#accuracy-scale");
    const asimResult = document.querySelector("#asim-result");
    const asimActive = candidate.architecture === "SRAM" && ["QS", "QR", "QS-QR"].includes(candidate.model);
    const asimFields = ["#asim-adc-boundary", "#asim-adc-only", "#asim-random", "#asim-nonlinearity"];
    title.textContent = asimActive ? "Network-level accuracy" : "Accuracy awareness";
    scale.hidden = asimActive;
    asimResult.hidden = !asimActive;
    if (asimActive) {
      marker.hidden = true;
      if (!Number.isInteger(candidate.inputBits) || !Number.isInteger(candidate.weightBits) || !Number.isInteger(candidate.adcBits)) {
        lead.textContent = "ASiM screening needs integer precisions and an ADC";
        detail.textContent = "Choose integer input and weight precision plus an integer ADC precision. ASiM decomposes quantized tensors into binary MAC cycles.";
        asimFields.forEach((selector) => { document.querySelector(selector).textContent = "Not evaluated"; });
        document.querySelector("#asim-evidence-class").textContent = "Outside ASiM input scope";
        document.querySelector("#asim-config-text").textContent = "Select integer input/weight precision and an ADC precision to generate an ASiM configuration.";
        return;
      }
      const screening = labAsimScreen(candidate);
      if (!screening) {
        lead.textContent = "ASiM tables are unavailable";
        detail.textContent = "The ASiM result tables did not load.";
        return;
      }
      lead.textContent = `${screening.bound ? "At most" : "About"} ${labFormat(screening.estimate, 1)}% top-1`;
      detail.textContent = `${screening.workload.label}, ${screening.encoding.label} encoding, scaled from the paper's ${labFormat(screening.paperBaseline, 2)}% digital baseline onto the entered ${labFormat(candidate.asimBaseline, 2)}%. ${screening.notes.join(" ")} This reads published ASiM results and is not an inference run.`;
      const noiseText = (source) => (source.sigma <= 0 ? "None" : `${labFormat(source.sigma, 3)}% Vpp (${labFormat(source.sigma / 100 * 255, 3)} LSB rms at 8 bits) → ${labFormat(source.scaled, 1)}% at an 8-bit ADC (${source.figure})`);
      document.querySelector("#asim-adc-boundary").textContent = `${screening.boundaryAdc}-bit boundary · ${labFormat(candidate.adcBits, 2)}-bit selected`;
      document.querySelector("#asim-adc-only").textContent = `${labFormat(screening.adcOnly, 1)}% (Figs. 6–7)`;
      document.querySelector("#asim-random").textContent = noiseText(screening.noise[0]);
      document.querySelector("#asim-nonlinearity").textContent = noiseText(screening.noise[1]);
      document.querySelector("#asim-evidence-class").textContent = screening.evidenceClass;
      document.querySelector("#asim-config-text").textContent = screening.config;
      return;
    }
    const presentSnr = candidate.hasSnr ? candidate.snr : modelResult.modeledSnr;
    if (presentSnr === null || !Number.isFinite(presentSnr)) {
      marker.hidden = true;
      if (modelResult.referenceSnr !== null && Number.isFinite(modelResult.referenceSnr)) {
        lead.textContent = "SNR not established at the selected node";
        detail.textContent = `The source model gives ${labFormat(modelResult.referenceSnr, 1)} dB at its ${modelResult.referenceNode} reference. A new technology node can increase nonidealities and reduce SNR, so this value is not transferred.`;
      } else {
        const gap = (candidate.adcBits || 0) - candidate.information;
        lead.textContent = "Accuracy not established";
        detail.textContent = candidate.adcBits === null ? "This compute model has no selected ADC. Add measured compute SNR or a validated architecture model before making an accuracy claim." : `ADC precision is ${Math.abs(gap).toFixed(1)} bits ${gap < 0 ? "below" : "above"} the selected pre-ADC information proxy. This flags quantization pressure but is not an accuracy prediction.`;
      }
      return;
    }
    marker.hidden = false;
    marker.style.left = `${Math.max(0, Math.min(100, presentSnr / 50 * 100))}%`;
    const sourceText = candidate.hasSnr ? "Measured" : "Model-estimated";
    if (presentSnr < 10) {
      lead.textContent = `${sourceText} SNR is below cited guidance`;
      detail.textContent = `${labFormat(presentSnr, 1)} dB is below the paper's task-dependent 10–40 dB context. Network accuracy requires direct evaluation.`;
    } else if (presentSnr <= 40) {
      lead.textContent = `${sourceText} SNR is within cited guidance`;
      detail.textContent = `${labFormat(presentSnr, 1)} dB falls inside the cited task-dependent range. This is context, not a network-accuracy guarantee.`;
    } else {
      lead.textContent = `${sourceText} SNR is above cited guidance`;
      detail.textContent = `${labFormat(presentSnr, 1)} dB exceeds the cited range, but mapping and target-network accuracy still require direct evaluation.`;
    }
  }

  function updateMetricCitations(candidate) {
    const modelCitation = document.querySelector("#model-citation");
    const accuracyCitation = document.querySelector("#accuracy-citation");
    if (candidate.architecture === "eNVM" && candidate.model === "IS") {
      modelCitation.innerHTML = '<span>Source</span><a href="https://doi.org/10.1109/JXCDC.2024.3381888">Roy &amp; Shanbhag, JXCDC 2024</a><a href="https://github.com/calmyor/eNVM-IMC-Modeling">Code</a>';
      accuracyCitation.innerHTML = '<span>Source</span><a href="https://doi.org/10.1109/JXCDC.2024.3381888">Resistive IMC SNDR model</a>';
      return;
    }
    if (candidate.architecture === "SRAM" && candidate.model === "QR") {
      modelCitation.innerHTML = '<span>Source</span><a href="https://doi.org/10.1109/TCSI.2025.3594230">Kavishwar &amp; Shanbhag, TCAS-I</a><a href="https://github.com/mihirvk2/tcas-mimo-imc-2025">Code</a>';
    } else if (candidate.architecture === "SRAM" && candidate.model === "DIMC") {
      modelCitation.innerHTML = '<span>Scope</span>Behavioral implementation reference; no paper-backed metric model';
    } else {
      modelCitation.innerHTML = '<span>Method</span><a href="https://doi.org/10.1109/CICC53496.2022.9772817">CICC benchmarking method</a>';
    }
    if (candidate.architecture === "SRAM" && ["QS", "QR", "QS-QR"].includes(candidate.model)) {
      accuracyCitation.innerHTML = '<span>Source</span><a href="https://arxiv.org/abs/2411.11022">Zhang et al., ASiM</a><a href="https://github.com/Keio-CSG/ASiM">Code</a>';
    } else if (candidate.model === "DIMC") {
      accuracyCitation.innerHTML = '<span>Scope</span>No network-accuracy source is attached';
    } else {
      accuracyCitation.innerHTML = '<span>Method</span><a href="https://doi.org/10.1109/OJSSCS.2022.3210152">OJ-SSCS accuracy framing</a>';
    }
  }

  function draw(candidate) {
    chart.replaceChildren(chart.querySelector("title"), chart.querySelector("desc"));
    const grid = svgElement("g");
    xTicks.forEach((tick) => {
      const x = xPosition(tick);
      grid.append(svgElement("line", { x1: x, y1: margin.top, x2: x, y2: margin.top + plotHeight, class: "grid-line" }));
      const label = svgElement("text", { x, y: margin.top + plotHeight + 28, "text-anchor": "middle" }); label.textContent = tick; grid.append(label);
    });
    yTicks.forEach((tick) => {
      const y = yPosition(tick);
      grid.append(svgElement("line", { x1: margin.left, y1: y, x2: margin.left + plotWidth, y2: y, class: "grid-line" }));
      const label = svgElement("text", { x: margin.left - 14, y: y + 5, "text-anchor": "end" }); label.textContent = tick.toLocaleString(); grid.append(label);
    });
    grid.append(svgElement("line", { x1: margin.left, y1: margin.top + plotHeight, x2: margin.left + plotWidth, y2: margin.top + plotHeight, class: "axis-line" }));
    grid.append(svgElement("line", { x1: margin.left, y1: margin.top, x2: margin.left, y2: margin.top + plotHeight, class: "axis-line" }));
    const xTitle = svgElement("text", { x: margin.left + plotWidth / 2, y: height - 18, "text-anchor": "middle", class: "axis-title" }); xTitle.textContent = "Compute density (1b-TOPS/mm²)"; grid.append(xTitle);
    const yTitle = svgElement("text", { x: 24, y: margin.top + plotHeight / 2, "text-anchor": "middle", class: "axis-title", transform: `rotate(-90 24 ${margin.top + plotHeight / 2})` }); yTitle.textContent = "Energy efficiency (1b-TOPS/W)"; grid.append(yTitle);
    chart.append(grid);
    const pointGroup = svgElement("g");
    plottedRows.forEach((row) => pointGroup.append(reportedPointLink(row, { cx: xPosition(row.density), cy: yPosition(row.efficiency), r: compact ? 5.2 : 4.2, fill: colors[row.Architecture] || "#687080" })));
    chart.append(pointGroup);
    const offChart = candidate.density < xDomain[0] || candidate.density > xDomain[1] || candidate.efficiency < yDomain[0] || candidate.efficiency > yDomain[1];
    const x = xPosition(Math.max(xDomain[0], Math.min(xDomain[1], candidate.density)));
    const y = yPosition(Math.max(yDomain[0], Math.min(yDomain[1], candidate.efficiency)));
    const size = compact ? 9 : 10;
    chart.append(svgElement("polygon", { points: `${x},${y-size} ${x+size},${y} ${x},${y+size} ${x-size},${y}`, fill: colors[candidate.architecture] || "#20242b", class: "candidate-point" }));
    const label = svgElement("text", { x: Math.min(x + 15, width - (offChart ? 190 : 118)), y: Math.max(y - 14, 24), class: "candidate-label" });
    label.textContent = `Your ${candidate.architecture}${offChart ? " (off chart)" : ""}`;
    chart.append(label);
  }

  function updateNeighbors(candidate) {
    const nearestByIndex = new Map();
    plottedRows.forEach((row) => {
      const distance = Math.hypot(Math.log10(row.density / candidate.density), Math.log10(row.efficiency / candidate.efficiency));
      const prior = nearestByIndex.get(row.Index);
      if (!prior || distance < prior.distance) nearestByIndex.set(row.Index, { ...row, distance });
    });
    const nearest = [...nearestByIndex.values()].sort((a, b) => a.distance - b.distance).slice(0, 3);
    document.querySelector("#neighbor-body").innerHTML = renderNearestPaperRows(nearest);
  }

  function update() {
    let candidate = candidateValues();
    const modelAvailable = Number.isFinite(candidate.modelEnergyPerOpFj);
    controls.energySource.querySelector('option[value="model"]').disabled = !modelAvailable;
    if (!modelAvailable && controls.energySource.value === "model") {
      controls.energySource.value = "entered";
      candidate = candidateValues();
    }
    controls.columnEnergy.disabled = candidate.usesModelEnergy;
    if (candidate.usesModelEnergy) {
      const logEnergy = Math.log10(candidate.columnEnergyFj);
      controls.columnEnergy.value = String(Math.max(Number(controls.columnEnergy.min), Math.min(Number(controls.columnEnergy.max), logEnergy)));
    }
    const opsPerColumn = 2 * candidate.dimension * candidate.inputBits * candidate.weightBits;
    document.querySelector("#lab-energy-source-help").textContent = !modelAvailable
      ? "No attached energy model covers this architecture and compute model."
      : candidate.usesModelEnergy
        ? `Using the model estimate: ${labFormat(candidate.modelEnergyPerOpFj, 3)} fJ per 1b op.`
        : `Model estimate available: ${labFormat(candidate.modelEnergyPerOpFj * opsPerColumn, 3)} fJ per column.`;

    outputs.tech.textContent = `${candidate.tech} nm`;
    outputs.coreLatency.textContent = `${labFormat(candidate.coreLatencyNs, 3)} ns`;
    outputs.area.textContent = `${labFormat(candidate.areaMm2, 4)} mm²`;
    outputs.columnEnergy.textContent = `${labFormat(candidate.columnEnergyFj, 3)} fJ`;
    outputs.information.textContent = `${labFormat(candidate.information, 1)} bits`;
    outputs.snr.textContent = `${labFormat(candidate.snr, 1)} dB`;
    outputs.vbl.textContent = `${labFormat(candidate.vbl, 4)} V`;
    outputs.cellCap.textContent = `${labFormat(candidate.cellCap, 3)} fF`;
    outputs.adcNoise.textContent = `${labFormat(candidate.adcNoiseMv, 2)} mV`;
    outputs.asimRandomNoise.innerHTML = `${labFormat(candidate.asimRandomNoise, 3)}% V<sub>pp</sub>`;
    outputs.asimNonlinearity.innerHTML = `${labFormat(candidate.asimNonlinearity, 3)}% V<sub>pp</sub>`;
    controls.snr.disabled = !candidate.hasSnr;
    if (candidate.envmDevice === "FeFET") controls.vbl.max = "0.7"; else controls.vbl.max = "-0.09";
    document.querySelector("#lab-derived-formula").textContent = `N₁b = 2 × ${candidate.adcColumns} × ${candidate.dimension} × ${labFormat(candidate.inputBits, 2)} × ${labFormat(candidate.weightBits, 2)} = ${labFormat(candidate.bitOps, 0)} per invocation · ${labFormat(candidate.energyPerOpFj, 3)} fJ per 1b op · ${labFormat(candidate.throughput, 3)} 1b-TOPS.`;

    const evidence = updateEvidence(candidate);
    const modelResult = validateModel(candidate);
    draw(candidate);
    document.querySelector("#placement-lead").textContent = `${candidate.architecture} · ${candidate.model || "unclassified"} · ${candidate.tech} nm`;
    document.querySelector("#derived-energy").textContent = `${labFormat(candidate.energyPerOpFj, 3)} fJ`;
    document.querySelector("#derived-throughput").textContent = `${labFormat(candidate.throughput, 3)} 1b-TOPS`;
    document.querySelector("#derived-efficiency").textContent = `${labFormat(candidate.efficiency, 1)} 1b-TOPS/W`;
    document.querySelector("#derived-density").textContent = `${labFormat(candidate.density, 2)} 1b-TOPS/mm²`;
    const comparison = evidence.pair.filter((row) => row.efficiency > 0 && row.density > 0);
    const efficiencyRange = labRange(comparison.map((row) => row.efficiency));
    const densityRange = labRange(comparison.map((row) => row.density));
    const inside = efficiencyRange && densityRange && candidate.efficiency >= efficiencyRange[0] && candidate.efficiency <= efficiencyRange[1] && candidate.density >= densityRange[0] && candidate.density <= densityRange[1];
    document.querySelector("#envelope-check").textContent = !efficiencyRange || !densityRange ? "No pair-specific envelope" : inside ? "Inside reported pair range" : "Outside reported pair range";
    document.querySelector("#model-lead").textContent = modelResult.lead;
    document.querySelector("#model-detail").textContent = modelResult.detail;
    document.querySelector("#model-metrics").innerHTML = modelResult.metrics.map(([label, value]) => `<li><span>${escapeHtml(label)}</span><strong>${escapeHtml(value)}</strong></li>`).join("");
    updateAccuracy(candidate, modelResult);
    updateMetricCitations(candidate);
    updateNeighbors(candidate);
  }

  let pendingUpdate = null;
  function scheduleUpdate() {
    clearTimeout(pendingUpdate);
    pendingUpdate = setTimeout(update, 90);
  }
  function applyAsimBaseline() {
    const tables = window.IMC_ASIM_DATA?.workloads?.[controls.asimWorkload.value];
    if (!tables) return;
    const baseline = controls.asimTraining.value === "qat" && tables.nonNatBaseline ? tables.nonNatBaseline : tables.baseline;
    controls.asimBaseline.value = String(baseline);
  }
  controls.architecture.addEventListener("change", () => { refreshModelOptions(); update(); });
  controls.model.addEventListener("change", () => { refreshModelScope(); update(); });
  controls.energySource.addEventListener("change", update);
  controls.envmDevice.addEventListener("change", () => {
    const maximum = controls.envmDevice.value === "FeFET" ? 0.7 : -0.09;
    controls.vbl.max = String(maximum);
    if (Number(controls.vbl.value) > maximum) controls.vbl.value = String(maximum);
    update();
  });
  controls.asimWorkload.addEventListener("change", () => { applyAsimBaseline(); update(); });
  controls.asimTraining.addEventListener("change", () => { applyAsimBaseline(); update(); });
  const handled = ["architecture", "model", "energySource", "envmDevice", "asimWorkload", "asimTraining"];
  Object.entries(controls).filter(([key]) => !handled.includes(key)).forEach(([, control]) => control.addEventListener("input", scheduleUpdate));
  applyAsimBaseline();
  refreshModelOptions();
  update();
}

initDetailedLab().catch((error) => {
  const lead = document.querySelector("#placement-lead");
  if (lead) lead.textContent = error.message;
});
