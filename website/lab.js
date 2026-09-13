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

// Compute-accuracy and energy models for charge- and current-domain SRAM IMC.
//   [G] S. K. Gonugondla, C. Sakr, H. Dbouk, N. R. Shanbhag, "Fundamental limits on energy-delay-accuracy of
//       in-memory architectures in inference applications," IEEE TCAD 2022 (arXiv:2012.13645):
//       Table II (65 nm parameters), Table III (architecture noise and energy), eqs. (10)-(11), (18), (21), (25).
//   [D] S. K. Roy, PhD dissertation, UIUC 2024: E_op1 = E_DP / (2·N·B_x·B_w) (2.3), SNDR_col (2.6), E_ADC (3.18).
// Operands follow [G]: unsigned inputs x in [0, 1) and signed weights w in [-1, 1), each built from
// independent equiprobable bits, which is what produces the (1 - 4^-B) factors in Table III.

const LAB_G_PARAMETERS = {
  kPrime: 220e-6, alpha: 1.8, sigmaVt: 23.8e-3, vt: 0.4,   // QS cell current, 65 nm
  kappa: 0.08, wlCox: 0.31, injection: 0.5,                // QR: Pelgrom coefficient (fF^0.5), switch W·L·Cox (fF), p
  boltzmann: 1.38e-23, temperature: 300
};

function labAdcEnergyFj(bits) {
  return 100 * bits + 1e-3 * 4 ** bits; // [D] (3.18): k1 = 100 fJ, k2 = 1 aJ
}

function labDb(ratio) {
  return ratio === Infinity ? Infinity : 10 * Math.log10(ratio);
}

function labParallel(...ratios) {
  const inverse = ratios.filter((ratio) => ratio !== null && ratio !== undefined).reduce((sum, ratio) => sum + (ratio === Infinity ? 0 : 1 / ratio), 0);
  return inverse === 0 ? Infinity : 1 / inverse;
}

function labOperands(inputBits, weightBits) {
  const xMean = (1 - 2 ** -inputBits) / 2;
  const xVariance = (1 - 4 ** -inputBits) / 12;
  const wMean = -(2 ** -weightBits);
  const wVariance = (1 - 4 ** -weightBits) / 3;
  const xSquare = xVariance + xMean ** 2;
  const wSquare = wVariance + wMean ** 2;
  return { xMean, xVariance, xSquare, wMean, wVariance, wSquare, productVariance: wSquare * xSquare - (wMean * xMean) ** 2 };
}

function labBinomial(n, p) {
  const values = [];
  const probs = [];
  let log = n * Math.log(1 - p);
  const odds = Math.log(p / (1 - p));
  for (let k = 0; k <= n; k += 1) {
    if (k > 0) log += Math.log((n - k + 1) / k) + odds;
    const probability = Math.exp(log);
    if (probability > 1e-14) {
      values.push(k);
      probs.push(probability);
    }
  }
  return { values, probs };
}

// Mid-rise uniform ADC: `levels` equal bins tile [low, high] and each input maps to its bin centre;
// inputs beyond the range fall into the end bins.
function labQuantize(value, low, step, levels) {
  return low + (Math.min(levels - 1, Math.max(0, Math.floor((value - low) / step))) + 0.5) * step;
}

function labDiscreteMse(distribution, low, high, levels) {
  const step = (high - low) / levels;
  let mse = 0;
  for (let i = 0; i < distribution.values.length; i += 1) {
    const error = distribution.values[i] - labQuantize(distribution.values[i], low, step, levels);
    mse += distribution.probs[i] * error * error;
  }
  return mse;
}

function labGaussianMse(mean, sd, low, high, levels) {
  const step = (high - low) / levels;
  const from = mean - 10 * sd;
  const span = 20 * sd;
  const samples = Math.min(200000, Math.max(4000, Math.ceil(span / (step / 8))));
  const width = span / samples;
  let mse = 0;
  let mass = 0;
  for (let i = 0; i < samples; i += 1) {
    const value = from + (i + 0.5) * width;
    const z = (value - mean) / sd;
    const density = Math.exp(-0.5 * z * z);
    const error = value - labQuantize(value, low, step, levels);
    mse += density * error * error;
    mass += density;
  }
  return mse / mass;
}

// Quantization MSE of a B-bit ADC with the clipping range that minimizes it (searched over the full input range
// and mean ± kσ, k = 1…8), in the units of the distribution. Returns zero error when every input level is resolved.
function labBestAdc(distribution, bits) {
  if (!(bits > 0)) return null;
  const levels = Math.max(2, Math.round(2 ** bits));
  let { mean, sd, min, max } = distribution;
  if (distribution.kind === "discrete") {
    let square = 0;
    mean = 0;
    min = Infinity;
    max = -Infinity;
    distribution.values.forEach((value, i) => {
      mean += distribution.probs[i] * value;
      square += distribution.probs[i] * value * value;
      min = Math.min(min, value);
      max = Math.max(max, value);
    });
    sd = Math.sqrt(Math.max(0, square - mean * mean));
  }
  if (distribution.lattice > 0 && Number.isFinite(min) && Number.isFinite(max) && Math.round((max - min) / distribution.lattice) + 1 <= levels) {
    return { mse: 0, lossless: true, low: min, high: max, clip: "every input level resolved" };
  }
  if (!(sd > 0)) return { mse: 0, lossless: true, low: mean, high: mean, clip: "constant input" };
  let best = null;
  const consider = (low, high, clip) => {
    if (!(high > low)) return;
    const mse = distribution.kind === "discrete" ? labDiscreteMse(distribution, low, high, levels) : labGaussianMse(mean, sd, low, high, levels);
    if (!best || mse < best.mse) best = { mse, low, high, clip, lossless: false };
  };
  if (Number.isFinite(min) && Number.isFinite(max)) consider(min, max, "full range");
  for (let k = 1; k <= 8.001; k += 0.25) consider(Math.max(min, mean - k * sd), Math.min(max, mean + k * sd), `±${k}σ`);
  return best;
}

// Pre-ADC value of one QR-Arch conversion: z = Σ_j x_j·b_j for one weight bit b_j, with B_x-bit unsigned x_j.
function labQrConversionDistribution(dimension, inputBits) {
  const codes = 2 ** inputBits;
  if (Number.isInteger(inputBits) && dimension * dimension * codes * codes / 2 <= 4e6) {
    const step = new Float64Array(codes).fill(0.5 / codes);
    step[0] += 0.5;
    let pmf = new Float64Array([1]);
    for (let row = 0; row < dimension; row += 1) {
      const next = new Float64Array(pmf.length + codes - 1);
      for (let i = 0; i < pmf.length; i += 1) {
        if (pmf[i] < 1e-300) continue;
        for (let code = 0; code < codes; code += 1) next[i + code] += pmf[i] * step[code];
      }
      pmf = next;
    }
    const values = [];
    const probs = [];
    pmf.forEach((probability, index) => {
      if (probability > 1e-14) {
        values.push(index / codes);
        probs.push(probability);
      }
    });
    return { kind: "discrete", values, probs, lattice: 1 / codes };
  }
  const x = labOperands(inputBits, 1);
  return {
    kind: "gaussian", mean: dimension * x.xMean / 2, sd: Math.sqrt(dimension * (x.xSquare / 2 - x.xMean ** 2 / 4)),
    min: 0, max: dimension * (1 - 2 ** -inputBits), lattice: Number.isInteger(inputBits) ? 1 / codes : 0
  };
}

function labQsCell(wordlineVoltage) {
  const p = LAB_G_PARAMETERS;
  const overdrive = wordlineVoltage - p.vt;
  if (!(overdrive > 0)) return null;
  return { current: p.kPrime * overdrive ** p.alpha, sigmaD: p.alpha * p.sigmaVt / overdrive }; // α-law current, (18)
}

// c: { arch: "QS-Arch" | "QR-Arch" | "CM" | "IS", dimension, inputBits, weightBits, adcBits, vdd,
//      wordlineVoltage (V), bitlineCapFf, pulsePs (LSB WL pulse), headroomV (ΔV_BL,max), unitCapFf (C_o), includeInjection }
function labChargeDomainModel(c) {
  const p = LAB_G_PARAMETERS;
  const N = Math.max(1, Math.round(c.dimension));
  const bx = c.inputBits;
  const bw = c.weightBits;
  const ops = labOperands(bx, bw);
  const signal = N * ops.productVariance;
  const inputQuantization = N * (4 ** -bx * ops.wVariance + 4 ** (1 - bw) * ops.xSquare) / 12; // σ²_qiy, Table III
  const notes = [];
  const result = { arch: c.arch, dimension: N, signal, inputSqnr: signal / inputQuantization, notes };
  const adcFj = c.adcBits > 0 ? labAdcEnergyFj(c.adcBits) : null;

  if (c.arch === "QS-Arch" || c.arch === "IS") {
    const cell = labQsCell(c.wordlineVoltage);
    if (!cell) return { ...result, error: `V_WL must exceed V_t = ${p.vt} V.` };
    const placeFactor = (1 - 4 ** -bw) * (1 - 4 ** -bx);
    const unitDischarge = cell.current * c.pulsePs * 1e-12 / (c.bitlineCapFf * 1e-15);
    const headroomUnits = c.arch === "IS" ? Infinity : c.headroomV / unitDischarge;
    const counts = labBinomial(N, 0.25);
    let clipping = 0;
    let meanUnits = 0;
    const clipped = counts.values.map((k, i) => {
      if (k > headroomUnits) clipping += counts.probs[i] * (k - headroomUnits) ** 2;
      meanUnits += counts.probs[i] * Math.min(k, headroomUnits);
      return Math.min(k, headroomUnits);
    });
    result.circuitNoise = N * cell.sigmaD ** 2 * placeFactor / 9;
    result.headroomNoise = 4 / 9 * placeFactor * clipping;
    result.analogSnr = signal / (result.circuitNoise + result.headroomNoise);
    result.sigmaD = cell.sigmaD;
    result.headroomUnits = headroomUnits;
    result.conversionsPerDp = bx * bw;
    const adc = labBestAdc({ kind: "discrete", values: clipped, probs: counts.probs, lattice: 1 }, c.adcBits);
    result.adc = adc;
    result.adcNoise = adc ? adc.mse * 4 / 9 * placeFactor : null;
    if (c.arch === "QS-Arch") {
      const bitlineFj = meanUnits * unitDischarge * c.vdd * c.bitlineCapFf; // (21): E[V_a]·V_dd·C_BL
      result.energy = adcFj === null ? null : { arrayFj: bitlineFj, adcFj, perDpFj: bx * bw * (bitlineFj + adcFj) };
      result.minAdcBits = (snrA) => Math.min((labDb(snrA) + 16.2) / 6, Math.log2(headroomUnits), Math.log2(N));
    } else {
      notes.push("No published current-summing SRAM model: SNR_a uses the QS-Arch cell-current mismatch term without headroom clipping, and no energy is estimated.");
      result.energy = null;
      result.minAdcBits = (snrA) => Math.min((labDb(snrA) + 16.2) / 6, Math.log2(N));
    }
  } else if (c.arch === "QR-Arch") {
    const capF = c.unitCapFf * 1e-15;
    const mismatch = ops.xSquare * p.kappa ** 2 / c.unitCapFf; // E[x²]·σ²_Co/C_o², σ_Co = κ√C_o
    const thermal = 2 * p.boltzmann * p.temperature / (capF * c.vdd ** 2);
    // Charge injection v_j = p·WLCox·(V_dd − V_t − V_j)/C_o (24) scales as 1/C_o, so its variance scales as 1/C_o².
    // Table III's caption prints E[x²]·WLCox/C_o; that form gives ~2 dB at 1 fF and +4.8 dB per 3× C_o, against
    // the paper's reported ~+8 dB, which the squared form reproduces.
    const injection = c.includeInjection ? ops.xSquare * (p.injection * p.wlCox / c.unitCapFf) ** 2 : 0;
    result.circuitNoise = 2 / 3 * (1 - 4 ** -bw) * N * (mismatch + thermal + injection);
    result.headroomNoise = 0;
    result.analogSnr = signal / result.circuitNoise;
    result.noiseShares = { mismatch, thermal, injection };
    result.conversionsPerDp = bw;
    const adc = labBestAdc(labQrConversionDistribution(N, bx), c.adcBits);
    result.adc = adc;
    result.adcNoise = adc ? adc.mse * 4 * (1 - 4 ** -bw) / 3 : null;
    const arrayFj = N * c.unitCapFf * c.vdd ** 2; // (25) + N·E_mult with V_j = x_j·V_dd: N·C_o·V_dd²
    result.energy = adcFj === null ? null : { arrayFj, adcFj, perDpFj: bw * (arrayFj + adcFj) };
    result.minAdcBits = (snrA) => Math.min((labDb(snrA) + 16.2) / 6, bx + Math.log2(N));
  } else if (c.arch === "CM") {
    const cell = labQsCell(c.wordlineVoltage);
    if (!cell) return { ...result, error: `V_WL must exceed V_t = ${p.vt} V.` };
    const unitDischarge = cell.current * c.pulsePs * 1e-12 / (c.bitlineCapFf * 1e-15);
    const headroomUnits = c.headroomV / unitDischarge;
    result.circuitNoise = 2 / 3 * N * ops.xSquare * (1 / 4 - 4 ** -bw) * cell.sigmaD ** 2;
    result.headroomNoise = N * ops.xSquare * ops.wVariance * headroomUnits ** -2 * 2 ** (2 * bw) * Math.max(0, 1 - 2 * headroomUnits * 2 ** -bw) ** 2 / 12;
    result.analogSnr = signal / (result.circuitNoise + result.headroomNoise);
    result.sigmaD = cell.sigmaD;
    result.headroomUnits = headroomUnits;
    result.conversionsPerDp = 1;
    const lattice = Number.isInteger(bx) && Number.isInteger(bw) ? 2 ** -bx * 2 ** (1 - bw) : 0;
    const adc = labBestAdc({ kind: "gaussian", mean: N * ops.wMean * ops.xMean, sd: Math.sqrt(signal), min: -N * (1 - 2 ** -bx), max: N * (1 - 2 ** (1 - bw)) * (1 - 2 ** -bx), lattice }, c.adcBits);
    result.adc = adc;
    result.adcNoise = adc ? adc.mse : null;
    const magnitudes = 2 ** Math.max(0, Math.round(bw) - 1);
    let meanUnits = 0;
    for (let code = -magnitudes; code < magnitudes; code += 1) meanUnits += Math.min(Math.abs(code), headroomUnits) / (2 * magnitudes);
    const bitlineFj = 2 * N * meanUnits * unitDischarge * c.vdd * c.bitlineCapFf;         // 2N·E_QS
    const aggregateFj = N * c.unitCapFf * c.vdd ** 2 * (1 - ops.xMean);                     // E_QR, V_j = x_j·V_dd
    const multiplierFj = ops.xMean * (1 - ops.wMean) * c.unitCapFf * c.vdd ** 2;             // E_mult
    result.energy = adcFj === null ? null : { arrayFj: bitlineFj + aggregateFj + multiplierFj, adcFj, perDpFj: bitlineFj + aggregateFj + multiplierFj + adcFj };
    result.minAdcBits = (snrA) => (labDb(snrA) + 16.2) / 6;
  } else {
    return { ...result, error: `Unknown architecture ${c.arch}` };
  }

  result.analogSnrDb = labDb(result.analogSnr);
  result.inputSqnrDb = labDb(result.inputSqnr);
  result.snrA = labParallel(result.analogSnr, result.inputSqnr);
  result.snrADb = labDb(result.snrA);
  result.minAdcBits = result.minAdcBits(result.snrA);
  if (result.adcNoise === null) {
    notes.push("No ADC precision is selected, so SQNR and SNR_T are not evaluated.");
    result.adcSqnrDb = null;
    result.totalSnrDb = null;
    result.totalFloatSnrDb = null;
  } else {
    result.adcSqnr = result.adcNoise === 0 ? Infinity : signal / result.adcNoise;
    result.adcSqnrDb = labDb(result.adcSqnr);
    result.totalSnrDb = labDb(labParallel(result.analogSnr, result.adcSqnr));                       // [D] (2.6), vs fixed point
    result.totalFloatSnrDb = labDb(labParallel(result.analogSnr, result.inputSqnr, result.adcSqnr)); // [G] (10)-(11)
  }
  if (result.energy) result.energy.perOpFj = result.energy.perDpFj / (2 * N * bx * bw); // [D] (2.3)
  return result;
}

const labAutoAnalogArch = { QS: "QS-Arch", QR: "QR-Arch", "QS-QR": "CM", IS: "IS" };
const labAnalogArchLabels = {
  "QS-Arch": "QS-Arch: binarized bit-serial dot products (B_x·B_w conversions each)",
  "QR-Arch": "QR-Arch: binary-weighted dot products (B_w conversions each)",
  CM: "Compute memory: one multi-bit conversion per dot product",
  IS: "Current summing, approximated with the QS-Arch cell-current mismatch term"
};

function labIsChargeDomain(candidate) {
  return ["SRAM", "eDRAM"].includes(candidate.architecture) && ["QS", "QR", "QS-QR", "IS"].includes(candidate.model);
}

function labResolvedAnalogArch(candidate) {
  return candidate.analogArch === "auto" ? labAutoAnalogArch[candidate.model] : candidate.analogArch;
}

let labChargeCache = { key: "", result: null };
function labChargeDomainCached(candidate) {
  const inputs = {
    arch: labResolvedAnalogArch(candidate), dimension: candidate.dimension, inputBits: candidate.inputBits, weightBits: candidate.weightBits,
    adcBits: candidate.adcBits, vdd: candidate.vdd, wordlineVoltage: candidate.wordlineVoltage, bitlineCapFf: candidate.bitlineCapFf,
    pulsePs: candidate.pulsePs, headroomV: candidate.headroomV, unitCapFf: candidate.cellCap, includeInjection: candidate.includeInjection
  };
  const key = JSON.stringify(inputs);
  if (labChargeCache.key !== key) labChargeCache = { key, result: labChargeDomainModel(inputs) };
  return labChargeCache.result;
}

// Published-model energy per dot product, normalized as E_op1 = E_DP / (2·N·B_x·B_w) ([D] eq. 2.3).
function labPublishedEnergy(candidate) {
  if (candidate.architecture === "eNVM" && candidate.model === "IS" && window.IMC_MODEL_DATA) {
    const model = labEnvmModel(candidate);
    if (!model) return null;
    const conversions = candidate.envmCore === "multibit" ? 1 : candidate.inputBits * candidate.weightBits;
    const perDpFj = conversions * (model.arrayFj + model.adcFj);
    return { perDpFj, perOpFj: perDpFj / (2 * model.modelDimension * candidate.inputBits * candidate.weightBits), arrayFj: model.arrayFj, adcFj: model.adcFj, conversions, model };
  }
  if (labIsChargeDomain(candidate)) {
    const model = labChargeDomainCached(candidate);
    if (!model || model.error || !model.energy) return null;
    return { perDpFj: model.energy.perDpFj, perOpFj: model.energy.perOpFj, arrayFj: model.energy.arrayFj, adcFj: model.energy.adcFj, conversions: model.conversionsPerDp, model };
  }
  return null;
}

function labFitEnergy(candidate) {
  const fit = labFitColumnEnergy(candidate);
  return fit ? { ...fit, perOpFj: fit.columnFj / (2 * candidate.dimension * candidate.inputBits * candidate.weightBits) } : null;
}

function labModelEnergyPerOp(candidate, source) {
  const estimate = source === "published" ? labPublishedEnergy(candidate) : source === "fit" ? labFitEnergy(candidate) : null;
  return estimate ? estimate.perOpFj : null;
}

function labComputeAccuracy(candidate) {
  if (labIsChargeDomain(candidate)) return { kind: "charge", model: labChargeDomainCached(candidate) };
  if (candidate.architecture === "eNVM" && candidate.model === "IS" && window.IMC_MODEL_DATA) return { kind: "resistive", model: labEnvmModel(candidate) };
  if (candidate.model === "DIMC" || candidate.architecture === "Digital") return { kind: "digital" };
  return { kind: "none" };
}

// Data-derived column energy from scripts/build-energy-fit.py (window.IMC_ENERGY_FIT).
// Returns null when no fitted family covers the selection or a required input is missing.
function labFitColumnEnergy(candidate) {
  const fit = window.IMC_ENERGY_FIT;
  if (!fit) return null;
  const family = candidate.model === "DIMC"
    ? "DIMC"
    : ["eNVM", "eFlash"].includes(candidate.architecture) && candidate.model === "IS"
      ? "resistive IS"
      : { IS: "SRAM IS", QR: "QR", QS: "QS", "QS-QR": "QS" }[candidate.model];
  if (!family || !(candidate.vdd > 0) || !(candidate.tech > 0)) return null;
  const scale = candidate.tech / fit.nodeReferenceNm;
  let columnFj;
  if (family === "DIMC") {
    columnFj = fit.dimc.aFj * (candidate.dimension * candidate.inputBits * candidate.weightBits) ** fit.dimc.opsExponent
      * candidate.vdd ** 2 * scale ** fit.dimc.nodeExponent;
  } else {
    if (!(candidate.adcBits > 0)) return null;
    columnFj = fit.analog.k1FjPerBit * candidate.adcBits * scale ** fit.analog.adcNodeExponent
      + fit.analog.arrayFjPerRowVolt2[family] * candidate.dimension * candidate.vdd ** 2 * scale ** fit.analog.arrayNodeExponent;
  }
  return { family, columnFj, validation: fit.validation.families[family] };
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
    tech: document.querySelector("#lab-tech"), vdd: document.querySelector("#lab-vdd"), dimension: document.querySelector("#lab-dimension"),
    inputBits: document.querySelector("#lab-input-bits"), weightBits: document.querySelector("#lab-weight-bits"),
    adcBits: document.querySelector("#lab-adc-bits"), adcColumns: document.querySelector("#lab-adc-columns"),
    coreLatency: document.querySelector("#lab-core-latency"), area: document.querySelector("#lab-area"),
    energySource: document.querySelector("#lab-energy-source"), columnEnergy: document.querySelector("#lab-column-energy"),
    information: document.querySelector("#lab-information"),
    hasSnr: document.querySelector("#lab-has-snr"), snr: document.querySelector("#lab-snr"),
    envmDevice: document.querySelector("#lab-envm-device"), vbl: document.querySelector("#lab-vbl"),
    cellCap: document.querySelector("#lab-cell-cap"), adcNoise: document.querySelector("#lab-adc-noise"),
    analogArch: document.querySelector("#lab-analog-arch"), wordlineVoltage: document.querySelector("#lab-vwl"),
    bitlineCap: document.querySelector("#lab-cbl"), pulse: document.querySelector("#lab-pulse"),
    headroom: document.querySelector("#lab-headroom"), injection: document.querySelector("#lab-injection"),
    envmCore: document.querySelector("#lab-envm-core"),
    asimWorkload: document.querySelector("#lab-asim-workload"), asimBaseline: document.querySelector("#lab-asim-baseline"),
    asimEncoding: document.querySelector("#lab-asim-encoding"), asimRandomNoise: document.querySelector("#lab-asim-random-noise"),
    asimNonlinearity: document.querySelector("#lab-asim-nonlinearity"), asimTraining: document.querySelector("#lab-asim-training")
  };
  const outputs = {
    tech: document.querySelector("#lab-tech-output"), coreLatency: document.querySelector("#lab-core-latency-output"),
    area: document.querySelector("#lab-area-output"), columnEnergy: document.querySelector("#lab-column-energy-output"),
    information: document.querySelector("#lab-information-output"), snr: document.querySelector("#lab-snr-output"),
    vbl: document.querySelector("#lab-vbl-output"), cellCap: document.querySelector("#lab-cell-cap-output"),
    adcNoise: document.querySelector("#lab-adc-noise-output"), wordlineVoltage: document.querySelector("#lab-vwl-output"),
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
      envmDevice: controls.envmDevice.value, vbl: 10 ** Number(controls.vbl.value), envmCore: controls.envmCore.value,
      cellCap: 10 ** Number(controls.cellCap.value), adcNoiseMv: Number(controls.adcNoise.value),
      vdd: positiveNumber(controls.vdd, 0.8), analogArch: controls.analogArch.value,
      wordlineVoltage: Number(controls.wordlineVoltage.value), bitlineCapFf: positiveNumber(controls.bitlineCap, 270),
      pulsePs: positiveNumber(controls.pulse, 100), headroomV: positiveNumber(controls.headroom, 0.8), includeInjection: controls.injection.checked,
      asimWorkload: controls.asimWorkload.value, asimBaseline: Number(controls.asimBaseline.value),
      asimEncoding: Number(controls.asimEncoding.value), asimRandomNoise: Number(controls.asimRandomNoise.value),
      asimNonlinearity: Number(controls.asimNonlinearity.value), asimTraining: controls.asimTraining.value
    };
    const opsPerColumn = 2 * candidate.dimension * candidate.inputBits * candidate.weightBits;
    candidate.bitOps = candidate.adcColumns * opsPerColumn;
    candidate.modelEnergyPerOpFj = labModelEnergyPerOp(candidate, candidate.energySource);
    candidate.usesModelEnergy = candidate.energySource !== "entered" && Number.isFinite(candidate.modelEnergyPerOpFj);
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
    document.querySelector("#charge-model-controls").hidden = !(["SRAM", "eDRAM"].includes(architecture) && ["QS", "QR", "QS-QR", "IS"].includes(model));
    document.querySelector("#asim-model-controls").hidden = !(architecture === "SRAM" && ["QS", "QR", "QS-QR"].includes(model));
  }

  function updateEvidence(candidate) {
    const pair = rows.filter((row) => row.Architecture === candidate.architecture && row.model === candidate.model);
    const node = pair.filter((row) => labEqual(row.tech, candidate.tech));
    const dimension = node.filter((row) => labEqual(row.dimension, candidate.dimension));
    const precision = dimension.filter((row) => labEqual(row.inputBits, candidate.inputBits) && labEqual(row.weightBits, candidate.weightBits));
    // Digital IMC has no ADC input here, so a DIMC row's B_ADC field (28 for #177) is not matched.
    const ignoresAdc = candidate.model === "DIMC" || candidate.architecture === "Digital";
    const sameAdc = (row) => ignoresAdc || (candidate.adcBits === null ? row.adcBits === null : labEqual(row.adcBits, candidate.adcBits));
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

  function ratioText(entered, estimate) {
    const ratio = entered / estimate;
    return `${labFormat(ratio >= 1 ? ratio : 1 / ratio, 2)}× ${ratio >= 1 ? "above" : "below"}`;
  }

  function updateEnergy(candidate) {
    const published = labPublishedEnergy(candidate);
    const fit = labFitEnergy(candidate);
    const placementSource = candidate.usesModelEnergy ? (candidate.energySource === "published" ? "published model" : "data fit") : "entered column energy";
    const metrics = [["Placement energy per 1b op", `${labFormat(candidate.energyPerOpFj, 3)} fJ (${placementSource})`]];
    const details = [];
    if (published) {
      metrics.push(
        ["Published model per 1b op", `${labFormat(published.perOpFj, 3)} fJ`],
        ["Published energy per dot product", `${labFormat(published.perDpFj, 1)} fJ = ${published.conversions} × (${labFormat(published.arrayFj, 1)} array + ${labFormat(published.adcFj, 1)} ADC)`]
      );
    }
    if (fit) metrics.push(["Data fit per 1b op", `${labFormat(fit.perOpFj, 3)} fJ (typical error ×${labFormat(fit.validation.typicalFactor, 1)}; 80% within ×${labFormat(fit.validation.p80Factor, 1)})`]);

    if (candidate.architecture === "eNVM" && candidate.model === "IS" && published) {
      const model = published.model;
      details.push(`${candidate.envmDevice} current-summing column at ${[
        model.dimensionClamped ? `N = ${labFormat(model.modelDimension, 0)} (the edge of the model's 16–5,000 range; ${candidate.dimension} selected)` : `N = ${candidate.dimension}`,
        `${labFormat(candidate.vbl, 4)} V bitline${model.vblClamped ? " (clamped to the published range)" : ""}`,
        `the source's fixed ${labFormat(model.periodNs, 0)} ns read period`,
        `a ${labFormat(model.adcBits, 2)}-bit ADC${model.assumedAdc ? " (assumed; none selected)" : ""}`,
        candidate.envmCore === "multibit" ? "one multi-bit conversion per dot product" : "binarized columns"
      ].join(", ")}. E_bl = V_DD·V_BL·E[G_eq]·T and E_ADC = 100 fJ·B + 1 aJ·4^B; drivers and bias are excluded.`);
    } else if (labIsChargeDomain(candidate)) {
      const arch = labResolvedAnalogArch(candidate);
      if (published) {
        details.push(`${labAnalogArchLabels[arch]}. Table III energy with Table II 65 nm parameters at ${labFormat(candidate.vdd, 2)} V, V_WL = ${labFormat(candidate.wordlineVoltage, 2)} V, C_BL = ${labFormat(candidate.bitlineCapFf, 0)} fF, C_o = ${labFormat(candidate.cellCap, 2)} fF, and E_ADC = 100 fJ·B + 1 aJ·4^B. Switch, DAC, and driver energy (E_su, E_misc) have no published values and are excluded.`);
      } else if (arch === "IS") {
        details.push("No published energy model exists for current-summing SRAM, so only the data fit is shown.");
      } else if (!(candidate.adcBits > 0)) {
        details.push("The published model needs an ADC precision.");
      }
      if (candidate.architecture === "eDRAM") details.push("eDRAM uses the SRAM charge-domain forms as an approximation.");
      if (candidate.tech !== 65) details.push(`Model parameters are for 65 nm; nothing is scaled to ${candidate.tech} nm.`);
    } else if (candidate.model === "DIMC") {
      details.push("No published DIMC energy model exists; the data fit scales with adder activity N·B_x·B_w.");
      metrics.push(
        ["Behavioral reference", "1024×1024 INT8 matmul with a forced 2,150 µs latency"],
        ["Reference throughput", `${labFormat(2 * 1024 ** 3 / 2150e-6 / 1e12, 2)} TOPS (not an energy or area model)`]
      );
    }
    if (fit) details.push(`The data fit covers the ${fit.family} family, validated leave-one-paper-out on ${fit.validation.rows} rows from ${fit.validation.papers} papers.`);

    const reference = published || fit;
    const lead = !reference
      ? "No energy model covers this selection"
      : candidate.usesModelEnergy
        ? `Placement uses the ${candidate.energySource === "published" ? "published model" : "data fit"}`
        : `Entered energy is ${ratioText(candidate.energyPerOpFj, reference.perOpFj)} the ${published ? "published model" : "data fit"}`;
    if (!reference) details.push("The lab can check reported coordinates and exact paper matches, but it cannot estimate energy for this architecture and compute model.");
    document.querySelector("#model-lead").textContent = lead;
    document.querySelector("#model-detail").textContent = details.join(" ");
    document.querySelector("#model-metrics").innerHTML = metrics.map(([label, value]) => `<li><span>${escapeHtml(label)}</span><strong>${escapeHtml(value)}</strong></li>`).join("");
  }

  function updateAccuracy(candidate) {
    const lead = document.querySelector("#accuracy-lead");
    const detail = document.querySelector("#accuracy-detail");
    const marker = document.querySelector("#accuracy-marker");
    const scale = document.querySelector("#accuracy-scale");
    const accuracy = labComputeAccuracy(candidate);
    const metrics = [];
    const details = [];
    let modelSnrDb = null;

    if (accuracy.kind === "charge") {
      const model = accuracy.model;
      const arch = labResolvedAnalogArch(candidate);
      if (model.error) {
        details.push(model.error);
      } else {
        modelSnrDb = model.totalSnrDb;
        const noiseTotal = model.circuitNoise + model.headroomNoise;
        metrics.push(["Analog SNR_a", `${labFormat(model.analogSnrDb, 1)} dB`]);
        if (model.headroomNoise > 0) metrics.push(["Analog noise split", `${labFormat(100 * model.circuitNoise / noiseTotal, 0)}% circuit · ${labFormat(100 * model.headroomNoise / noiseTotal, 0)}% headroom clipping`]);
        metrics.push(["ADC SQNR", model.adcSqnrDb === null ? "No ADC selected" : model.adc.lossless ? "Lossless: every level resolved" : `${labFormat(model.adcSqnrDb, 1)} dB (optimal clip ${model.adc.clip})`]);
        metrics.push(["SNR_T vs fixed point", model.totalSnrDb === null ? "Needs an ADC" : `${labFormat(model.totalSnrDb, 1)} dB`]);
        if (model.totalFloatSnrDb !== null) metrics.push(["With operand quantization", `${labFormat(model.totalFloatSnrDb, 1)} dB (SQNR_qiy ${labFormat(model.inputSqnrDb, 1)} dB)`]);
        const minimumAdc = Math.max(1, model.minAdcBits);
        metrics.push(["Minimum ADC precision (MPC)", `${labFormat(minimumAdc, 1)} bits${candidate.adcBits ? (candidate.adcBits + 1e-9 >= minimumAdc ? " · selected ADC meets it" : " · selected ADC is below it") : ""}`]);
        if (Number.isFinite(model.headroomUnits)) metrics.push(["Bitline headroom", `${labFormat(model.headroomUnits, 1)} unit discharges · σ_I/I ${labFormat(100 * model.sigmaD, 1)}%`]);
        details.push(`${labAnalogArchLabels[arch]}. Table III noise with Table II 65 nm parameters and uniform, independent operand bits; the ADC SQNR is exact for the digitized signal at its SQNR-optimal clipping range.`);
        details.push(...model.notes);
        if (candidate.architecture === "eDRAM") details.push("eDRAM uses the SRAM charge-domain forms as an approximation.");
        if (candidate.tech !== 65) details.push(`No scaling to ${candidate.tech} nm is applied.`);
      }
      if (candidate.architecture === "SRAM" && candidate.model === "QR") {
        const monteCarlo = labSramQrSnrCached(candidate);
        const text = (result) => (result.errorFree ? "no errors" : `${labFormat(result.db, 1)} dB`);
        if (monteCarlo) metrics.push(["28 nm Monte Carlo cross-check", `${text(monteCarlo.fullRange)} full-range${monteCarlo.tunedWindow ? ` · ${text(monteCarlo.tunedWindow)} tuned` : ""} (${monteCarlo.samples} samples${monteCarlo.dimension < monteCarlo.requestedDimension ? ", N capped at 1,024" : ""})`]);
      }
    } else if (accuracy.kind === "resistive") {
      const model = accuracy.model;
      modelSnrDb = model.sndr.db;
      metrics.push(
        ["SNDR, 6-bit surface", `${labFormat(model.sndrSixBit, 1)} dB`],
        [`SNDR at ${labFormat(model.adcBits, 2)}-bit ADC`, model.sndr.db === null ? `Not modeled: ${model.sndr.reason}` : `${labFormat(model.sndr.db, 1)} dB${model.sndr.method === "adc-sweep" ? " (adjusted with the N = 64 ADC sweep)" : ""}`]
      );
      details.push(`${candidate.envmDevice} current-summing behavioral model: conductance variation, wire parasitics, current-mirror mismatch, and ADC noise at the SNDR-optimal clipping. Calibrated against a 22 nm MRAM prototype, with no node scaling.`);
    } else if (accuracy.kind === "digital") {
      modelSnrDb = Infinity;
      metrics.push(["Analog noise", "None: digital computation"], ["SNR_T vs fixed point", "Exact"]);
      details.push("Digital IMC computes the fixed-point dot product exactly, so accuracy is set by operand precision alone.");
    }
    document.querySelector("#accuracy-metrics").innerHTML = metrics.map(([label, value]) => `<li><span>${escapeHtml(label)}</span><strong>${escapeHtml(value)}</strong></li>`).join("");

    const presentSnr = candidate.hasSnr ? candidate.snr : modelSnrDb;
    const sourceText = candidate.hasSnr ? "Measured" : "Model";
    if (presentSnr === Infinity) {
      marker.hidden = true;
      scale.hidden = true;
      lead.textContent = "Exact fixed-point computation";
    } else if (presentSnr === null || !Number.isFinite(presentSnr)) {
      marker.hidden = true;
      scale.hidden = false;
      if (accuracy.kind === "charge" && !accuracy.model.error) {
        lead.textContent = "Select an ADC to evaluate SNR_T";
      } else if (accuracy.kind === "resistive") {
        lead.textContent = "SNDR not modeled for this ADC setting";
      } else {
        const gap = (candidate.adcBits || 0) - candidate.information;
        lead.textContent = "Accuracy not established";
        details.push(candidate.adcBits === null ? "No model covers this selection. Add a measured compute SNR before making an accuracy claim." : `No model covers this selection. ADC precision is ${Math.abs(gap).toFixed(1)} bits ${gap < 0 ? "below" : "above"} the pre-ADC information proxy, which flags quantization pressure but does not predict accuracy.`);
      }
    } else {
      scale.hidden = false;
      marker.hidden = false;
      marker.style.left = `${Math.max(0, Math.min(100, presentSnr / 50 * 100))}%`;
      lead.textContent = `${sourceText} SNR_T ${labFormat(presentSnr, 1)} dB`;
      details.push(presentSnr < 10
        ? "That is below the 10–40 dB range that fixed-point DNNs need to stay within 1% of floating point; network accuracy needs direct evaluation."
        : presentSnr <= 40
          ? "That falls inside the task-dependent 10–40 dB range reported for fixed-point DNNs within 1% of floating point; it is context, not a network-accuracy guarantee."
          : "That exceeds the 10–40 dB range reported for fixed-point DNNs; mapping and network accuracy still need direct evaluation.");
    }
    detail.textContent = details.join(" ");
    updateAsim(candidate);
  }

  function updateAsim(candidate) {
    const asimCard = document.querySelector("#asim-card");
    const asimActive = candidate.architecture === "SRAM" && ["QS", "QR", "QS-QR"].includes(candidate.model);
    asimCard.hidden = !asimActive;
    document.querySelector(".lab-results").classList.toggle("four-up", asimActive);
    if (!asimActive) return;
    const lead = document.querySelector("#asim-lead");
    const detail = document.querySelector("#asim-detail");
    const workload = document.querySelector("#asim-workload");
    workload.textContent = "";
    const asimFields = ["#asim-adc-boundary", "#asim-adc-only", "#asim-random", "#asim-nonlinearity"];
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
    workload.textContent = screening.workload.label;
    detail.textContent = `${screening.workload.label}, ${screening.encoding.label} encoding, scaled from the paper's ${labFormat(screening.paperBaseline, 2)}% digital baseline onto the entered ${labFormat(candidate.asimBaseline, 2)}%. ${screening.notes.join(" ")} This reads published ASiM results and is not an inference run.`;
    const noiseText = (source) => (source.sigma <= 0 ? "None" : `${labFormat(source.sigma, 3)}% Vpp (${labFormat(source.sigma / 100 * 255, 3)} LSB rms at 8 bits) → ${labFormat(source.scaled, 1)}% at an 8-bit ADC (${source.figure})`);
    document.querySelector("#asim-adc-boundary").textContent = `${screening.boundaryAdc}-bit boundary · ${labFormat(candidate.adcBits, 2)}-bit selected`;
    document.querySelector("#asim-adc-only").textContent = `${labFormat(screening.adcOnly, 1)}% (Figs. 6–7)`;
    document.querySelector("#asim-random").textContent = noiseText(screening.noise[0]);
    document.querySelector("#asim-nonlinearity").textContent = noiseText(screening.noise[1]);
    document.querySelector("#asim-evidence-class").textContent = screening.evidenceClass;
    document.querySelector("#asim-config-text").textContent = screening.config;
  }

  function updateMetricCitations(candidate) {
    const modelCitation = document.querySelector("#model-citation");
    const accuracyCitation = document.querySelector("#accuracy-citation");
    const gonugondla = '<a href="https://arxiv.org/abs/2012.13645">Gonugondla et al., TCAD 2022</a>';
    const fitLink = '<a href="data/Benchmarking_Data.csv">Data fit: Benchmarking_Data.csv</a>';
    if (candidate.architecture === "eNVM" && candidate.model === "IS") {
      modelCitation.innerHTML = `<span>Source</span><a href="https://doi.org/10.1109/JXCDC.2024.3381888">Roy &amp; Shanbhag, JXCDC 2024</a><a href="https://github.com/calmyor/eNVM-IMC-Modeling">Code</a>${fitLink}`;
      accuracyCitation.innerHTML = '<span>Source</span><a href="https://doi.org/10.1109/JXCDC.2024.3381888">Resistive IMC SNDR model</a>';
    } else if (labIsChargeDomain(candidate)) {
      modelCitation.innerHTML = `<span>Source</span>${gonugondla}${fitLink}`;
      accuracyCitation.innerHTML = `<span>Source</span>${gonugondla}${candidate.architecture === "SRAM" && candidate.model === "QR" ? '<a href="https://github.com/mihirvk2/tcas-mimo-imc-2025">28 nm QR code</a>' : ""}`;
    } else if (candidate.model === "DIMC" || candidate.architecture === "Digital") {
      modelCitation.innerHTML = `<span>Source</span>${fitLink}`;
      accuracyCitation.innerHTML = "<span>Scope</span>Digital computation has no analog noise";
    } else {
      modelCitation.innerHTML = `<span>Method</span><a href="https://doi.org/10.1109/CICC53496.2022.9772817">CICC benchmarking method</a>${labFitColumnEnergy(candidate) ? fitLink : ""}`;
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

  const labTrainingLabels = {
    nat: "Noise-aware or in-situ training",
    qat: "Quantization-aware training only",
    unstated: "Training not stated"
  };

  function networkAccuracyRows(indices) {
    const papers = window.IMC_NETWORK_ACCURACY.papers;
    return indices.map((index) => {
      const paper = papers[index];
      const row = rows.find((item) => item.Index === index);
      const link = window.BENCHMARK_PAPER_LINKS?.[index];
      const title = link
        ? `<a class="neighbor-paper-link" href="${escapeHtml(link.url)}" target="_blank" rel="noopener">${escapeHtml(row["Paper Title"])} ↗</a>`
        : escapeHtml(row["Paper Title"]);
      const results = paper.entries.length
        ? paper.entries.map((entry) => `<div>${escapeHtml(entry.task)}: <strong>${escapeHtml(entry.chipText || `${labFormat(entry.chip, 2)}%`)}</strong>${entry.baseline === null ? "" : ` vs ${labFormat(entry.baseline, 2)}% ${escapeHtml(entry.baselineType)}`} · ${escapeHtml(entry.evidence)}</div>`).join("")
        : `<div>${escapeHtml(paper.summary)}</div>`;
      const sources = paper.sources.map((source) => `<a class="neighbor-source-link" href="${escapeHtml(source.url)}" target="_blank" rel="noopener">${escapeHtml(source.label)} ↗</a>`).join(" · ");
      return `<tr><td><span class="index-chip">#${escapeHtml(index)}</span><div class="paper-meta">${escapeHtml(row.Year)} · ${escapeHtml(row["Journal/Conference"])}</div></td>`
        + `<td>${title}${results}<div class="paper-meta">Confidence: ${escapeHtml(paper.confidence)} · ${sources}</div></td>`
        + `<td><strong>${escapeHtml(labTrainingLabels[paper.training])}</strong><div>${escapeHtml(paper.trainingNote)}</div>${paper.compensation ? `<div class="paper-meta">${escapeHtml(paper.compensation)}</div>` : ""}</td>`
        + `<td>${escapeHtml(paper.bankMetric || "—")}</td></tr>`;
    }).join("");
  }

  function updateNetworkAccuracy(candidate, evidence) {
    const lead = document.querySelector("#network-lead");
    const detail = document.querySelector("#network-detail");
    const papers = window.IMC_NETWORK_ACCURACY?.papers;
    if (!papers) {
      lead.textContent = "Reported accuracy data did not load";
      return;
    }
    const familyIndices = [...new Set(rows.filter((row) => row.Architecture === candidate.architecture && row.model === candidate.model).map((row) => row.Index))];
    const exact = new Set(labUniquePapers(evidence.exact));
    const withData = familyIndices.filter((index) => papers[index]).sort((a, b) => (exact.has(b) - exact.has(a)) || Number(a) - Number(b));
    const trained = withData.filter((index) => papers[index].training === "nat");
    const family = `${candidate.architecture} ${candidate.model || "unclassified"}`;
    lead.textContent = withData.length
      ? `${withData.length} of ${familyIndices.length} ${family} papers report network accuracy`
      : `No ${family} paper reports network accuracy`;

    const sentences = [];
    if (trained.length) sentences.push(`${trained.map((index) => `#${index}`).join(", ")} used noise-aware or in-situ training.`);
    const exactWithout = [...exact].filter((index) => !papers[index]);
    if (exactWithout.length) sentences.push(`The matching paper${exactWithout.length === 1 ? "" : "s"} ${exactWithout.map((index) => `#${index}`).join(", ")} report${exactWithout.length === 1 ? "s" : ""} no network accuracy.`);
    const accuracy = labComputeAccuracy(candidate);
    const snr = candidate.hasSnr
      ? candidate.snr
      : accuracy.kind === "charge" ? accuracy.model.totalSnrDb
        : accuracy.kind === "resistive" ? accuracy.model.sndr.db
          : accuracy.kind === "digital" ? Infinity : null;
    if (snr === Infinity) {
      sentences.push("Digital computation adds no analog error, so reported accuracy reflects quantization alone.");
    } else if (Number.isFinite(snr) && snr < 20) {
      sentences.push(`At ${labFormat(snr, 1)} dB SNR_T, the benchmark's evidence is that near-baseline accuracy at low bank-level SNR came with noise-aware or in-situ training (#38, #13, #96) or was shown only for binary networks on MNIST. Binary networks lost 1.4–5.1 points on CIFAR-10 (#6, #17, #40), and the MRAM chip in Roy's dissertation, without noise-aware training, reached 74.8–82.0% against a 91.1% baseline at 5–11 dB SNDR.`);
    } else if (Number.isFinite(snr)) {
      sentences.push(`At ${labFormat(snr, 1)} dB SNR_T, the closest evidence is the high-SNR charge-domain chips that stayed within 0.6 points of quantized baselines without noise-aware training (#5, #25, #35; #25 and #35 measured 0.68 and 0.98 LSB rms column noise).`);
    }
    sentences.push("Each paper reports against its own baseline, often quantized software, and some results are simulated from measured chip statistics.");
    detail.textContent = sentences.join(" ");
    document.querySelector("#network-body").innerHTML = withData.length
      ? networkAccuracyRows(withData)
      : '<tr><td colspan="4">No reported network accuracy was found for this architecture and compute model.</td></tr>';
    const others = Object.keys(papers).filter((index) => !familyIndices.includes(index)).sort((a, b) => Number(a) - Number(b));
    document.querySelector("#network-other-summary").textContent = `${others.length} other analog papers with reported accuracy`;
    document.querySelector("#network-other-body").innerHTML = networkAccuracyRows(others);
  }

  function update() {
    let candidate = candidateValues();
    const available = { published: Number.isFinite(labModelEnergyPerOp(candidate, "published")), fit: Number.isFinite(labModelEnergyPerOp(candidate, "fit")) };
    ["published", "fit"].forEach((source) => { controls.energySource.querySelector(`option[value="${source}"]`).disabled = !available[source]; });
    if (controls.energySource.value !== "entered" && !available[controls.energySource.value]) {
      controls.energySource.value = "entered";
      candidate = candidateValues();
    }
    controls.columnEnergy.disabled = candidate.usesModelEnergy;
    if (candidate.usesModelEnergy) {
      const logEnergy = Math.log10(candidate.columnEnergyFj);
      controls.columnEnergy.value = String(Math.max(Number(controls.columnEnergy.min), Math.min(Number(controls.columnEnergy.max), logEnergy)));
    }
    const opsPerColumn = 2 * candidate.dimension * candidate.inputBits * candidate.weightBits;
    const columnEstimates = ["published", "fit"].filter((source) => available[source]).map((source) => `${source === "published" ? "published model" : "data fit"} ${labFormat(labModelEnergyPerOp(candidate, source) * opsPerColumn, 3)} fJ`);
    document.querySelector("#lab-energy-source-help").textContent = !columnEstimates.length
      ? "No energy model covers this architecture and compute model."
      : candidate.usesModelEnergy
        ? `Using the ${candidate.energySource === "published" ? "published model" : "data fit"}: ${labFormat(candidate.energyPerOpFj, 3)} fJ per 1b op.`
        : `Column estimates: ${columnEstimates.join(" · ")}.`;

    outputs.tech.textContent = `${candidate.tech} nm`;
    outputs.coreLatency.textContent = `${labFormat(candidate.coreLatencyNs, 3)} ns`;
    outputs.area.textContent = `${labFormat(candidate.areaMm2, 4)} mm²`;
    outputs.columnEnergy.textContent = `${labFormat(candidate.columnEnergyFj, 3)} fJ`;
    outputs.information.textContent = `${labFormat(candidate.information, 1)} bits`;
    outputs.snr.textContent = `${labFormat(candidate.snr, 1)} dB`;
    outputs.vbl.textContent = `${labFormat(candidate.vbl, 4)} V`;
    outputs.cellCap.textContent = `${labFormat(candidate.cellCap, 3)} fF`;
    outputs.adcNoise.textContent = `${labFormat(candidate.adcNoiseMv, 2)} mV`;
    outputs.wordlineVoltage.textContent = `${labFormat(candidate.wordlineVoltage, 2)} V`;
    outputs.asimRandomNoise.innerHTML = `${labFormat(candidate.asimRandomNoise, 3)}% V<sub>pp</sub>`;
    outputs.asimNonlinearity.innerHTML = `${labFormat(candidate.asimNonlinearity, 3)}% V<sub>pp</sub>`;
    controls.snr.disabled = !candidate.hasSnr;
    if (candidate.envmDevice === "FeFET") controls.vbl.max = "0.7"; else controls.vbl.max = "-0.09";
    document.querySelector("#lab-derived-formula").textContent = `N₁b = 2 × ${candidate.adcColumns} × ${candidate.dimension} × ${labFormat(candidate.inputBits, 2)} × ${labFormat(candidate.weightBits, 2)} = ${labFormat(candidate.bitOps, 0)} per invocation · ${labFormat(candidate.energyPerOpFj, 3)} fJ per 1b op · ${labFormat(candidate.throughput, 3)} 1b-TOPS.`;

    const evidence = updateEvidence(candidate);
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
    updateEnergy(candidate);
    updateAccuracy(candidate);
    updateMetricCitations(candidate);
    updateNeighbors(candidate);
    updateNetworkAccuracy(candidate, evidence);
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
