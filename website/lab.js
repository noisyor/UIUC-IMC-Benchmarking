const labModelLabels = {
  QS: "Charge sharing (QS)",
  QR: "Charge redistribution (QR)",
  DIMC: "Digital IMC (DIMC)",
  IS: "Current summing (IS)",
  "QS-QR": "Hybrid charge sharing / redistribution (QS–QR)"
};

const labModelOrder = ["QR", "QS", "DIMC", "IS", "QS-QR"];

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

function labBracket(values, selected, logarithmic = false) {
  const transform = logarithmic ? Math.log : (value) => value;
  const target = Math.max(values[0], Math.min(values[values.length - 1], selected));
  let upper = values.findIndex((value) => value >= target);
  if (upper === -1) return { lower: values.length - 1, upper: values.length - 1, mix: 0, clamped: true };
  if (upper === 0) return { lower: 0, upper: 0, mix: 0, clamped: target !== selected };
  const lower = upper - 1;
  const mix = (transform(target) - transform(values[lower])) / (transform(values[upper]) - transform(values[lower]));
  return { lower, upper, mix, clamped: target !== selected };
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

function labSignedBits(value, bits) {
  const modulus = 2 ** bits;
  let encoded = value < 0 ? modulus + value : value;
  return Array.from({ length: bits }, (_, index) => (encoded >> (bits - index - 1)) & 1);
}

function labPots(values) {
  return values.reduce((sum, value, index) => sum + value * (index === 0 ? -(2 ** (values.length - 1)) : 2 ** (values.length - index - 1)), 0);
}

function labVariance(values) {
  const mean = values.reduce((sum, value) => sum + value, 0) / values.length;
  return values.reduce((sum, value) => sum + (value - mean) ** 2, 0) / values.length;
}

function labSramQrSnr(candidate) {
  if (!Number.isInteger(candidate.inputBits) || !Number.isInteger(candidate.weightBits) || !Number.isInteger(candidate.adcBits)) return null;
  const dimension = Math.max(1, Math.min(1024, Math.round(candidate.dimension)));
  const sampleCount = Math.max(24, Math.min(96, Math.floor(12000 / dimension)));
  const random = labSeededRandom(280000 + dimension * 31 + candidate.inputBits * 7 + candidate.weightBits * 13 + candidate.adcBits * 17);
  const cMean = candidate.cellCap;
  const cSigma = 2.1 * 10 ** -2.5 * Math.sqrt(cMean);
  const vdd = 0.9;
  const cPar = cMean * dimension * 0.3 + 2.04278;
  const delta = cMean * vdd / (dimension * cMean + cPar);
  const levels = 2 ** candidate.adcBits;
  if (levels <= 2) return null;
  const t1 = 0.5 * delta * dimension / levels;
  const tM = (levels - 1.5) * delta * dimension / levels;
  const thresholdCount = levels - 1;
  const vLsb = (tM - t1) / (thresholdCount - 1);
  const thresholds = Array.from({ length: thresholdCount }, (_, index) => t1 + (tM - t1) * index / (thresholdCount - 1));
  const representations = [
    (t1 - vLsb / 2) / delta,
    ...thresholds.map((threshold) => (threshold + vLsb / 2) / delta)
  ];
  const capacitance = Array.from({ length: dimension }, () => Array.from({ length: candidate.weightBits }, () => Math.max(1e-6, cMean + cSigma * labNormal(random))));
  const capTotals = Array.from({ length: candidate.weightBits }, (_, bit) => cPar + capacitance.reduce((sum, row) => sum + row[bit], 0));
  const ideal = [];
  const modeled = [];

  function convert(voltage) {
    const noisy = voltage + candidate.adcNoiseMv * 1e-3 * labNormal(random);
    if (noisy < t1) return representations[0];
    if (noisy >= thresholds[thresholds.length - 1]) return representations[representations.length - 1];
    for (let index = 0; index < thresholds.length; index += 1) {
      if (noisy < thresholds[index]) return representations[index];
    }
    return representations[representations.length - 1];
  }

  for (let sample = 0; sample < sampleCount; sample += 1) {
    const weights = Array.from({ length: dimension }, () => Math.floor(random() * 2 ** candidate.weightBits) - 2 ** (candidate.weightBits - 1));
    const inputs = Array.from({ length: dimension }, () => Math.floor(random() * 2 ** candidate.inputBits) - 2 ** (candidate.inputBits - 1));
    const weightBits = weights.map((value) => labSignedBits(value, candidate.weightBits));
    const inputBits = inputs.map((value) => labSignedBits(value, candidate.inputBits));
    ideal.push(weights.reduce((sum, weight, index) => sum + weight * inputs[index], 0));
    const inputColumns = [];
    for (let inputBit = 0; inputBit < candidate.inputBits; inputBit += 1) {
      const weightColumns = [];
      for (let weightBit = 0; weightBit < candidate.weightBits; weightBit += 1) {
        let activeCap = 0;
        for (let row = 0; row < dimension; row += 1) activeCap += weightBits[row][weightBit] * inputBits[row][inputBit] * capacitance[row][weightBit];
        weightColumns.push(convert(activeCap * vdd / capTotals[weightBit]));
      }
      inputColumns.push(labPots(weightColumns));
    }
    modeled.push(labPots(inputColumns));
  }
  const errors = ideal.map((value, index) => value - modeled[index]);
  const errorVariance = labVariance(errors);
  return errorVariance > 0 ? 10 * Math.log10(labVariance(ideal) / errorVariance) : 100;
}

const labAsimWorkloads = {
  "resnet18-cifar10": { label: "ResNet-18 · CIFAR-10", task: "CIFAR-10", model: "ResNet-18", chance: 10, adcOffset: { 0: -2, 2: -1, 4: -1 }, noiseLimit: { 0: 0.4, 2: 0.2, 4: 0.1 }, sensitivity: "moderate" },
  "resnet18-imagenet": { label: "ResNet-18 · ImageNet", task: "ImageNet", model: "ResNet-18", chance: 0.1, adcOffset: { 0: 0, 2: 1, 4: 2 }, noiseLimit: { 0: 0.125, 2: 0.09, 4: 0.05 }, sensitivity: "high" },
  "vitb32-cifar10": { label: "ViT-B-32 · CIFAR-10", task: "CIFAR-10", model: "ViT-B-32", chance: 10, adcOffset: { 0: 0, 2: 1, 4: 2 }, noiseLimit: { 0: 0.125, 2: 0.09, 4: 0.05 }, sensitivity: "high" },
  "vitb32-imagenet": { label: "ViT-B-32 · ImageNet", task: "ImageNet", model: "ViT-B-32", chance: 0.1, adcOffset: { 0: 0, 2: 2, 4: 3 }, noiseLimit: { 0: 0.125, 2: 0.075, 4: 0.05 }, sensitivity: "very high" }
};

function labAsimScreen(candidate) {
  const workload = labAsimWorkloads[candidate.asimWorkload];
  const encoding = candidate.asimEncoding;
  const baseAdc = Math.ceil(Math.log2(Math.max(2, candidate.dimension + 1)));
  const fullRangeAdc = baseAdc + encoding;
  const requiredAdc = Math.max(1, baseAdc + workload.adcOffset[encoding]);
  const adcDeficit = Math.max(0, requiredAdc - (candidate.adcBits || 0));
  const combinedNoisePct = Math.hypot(candidate.asimRandomNoise, candidate.asimNonlinearity);
  const noiseLsb = candidate.adcBits === null ? Infinity : combinedNoisePct / 100 * (2 ** candidate.adcBits - 1);
  const noiseRatio = noiseLsb / workload.noiseLimit[encoding];

  let adcLoss;
  if (adcDeficit === 0) adcLoss = workload.sensitivity === "moderate" ? [0, 1.5] : [0, 3];
  else if (workload.sensitivity === "moderate") adcLoss = adcDeficit === 1 ? [1, 5] : [5 + 6 * (adcDeficit - 2), 16 + 10 * (adcDeficit - 2)];
  else if (workload.model === "ResNet-18") adcLoss = adcDeficit === 1 ? [7, 13] : [18 + 10 * (adcDeficit - 2), 40 + 15 * (adcDeficit - 2)];
  else adcLoss = adcDeficit === 1 ? [15, 40] : [35 + 15 * (adcDeficit - 2), 90];

  let noiseLoss;
  if (noiseRatio <= 0.5) noiseLoss = [0, 1];
  else if (noiseRatio <= 1) noiseLoss = [1, workload.sensitivity === "moderate" ? 3 : 5];
  else if (noiseRatio <= 2) noiseLoss = workload.sensitivity === "moderate" ? [3, 12] : [8, 25];
  else if (noiseRatio <= 4) noiseLoss = workload.sensitivity === "moderate" ? [10, 30] : [20, 60];
  else noiseLoss = workload.sensitivity === "moderate" ? [25, 70] : [45, 100];

  const minimumLoss = adcLoss[0] + noiseLoss[0];
  const maximumLoss = adcLoss[1] + noiseLoss[1];
  const low = Math.max(workload.chance, candidate.asimBaseline - maximumLoss);
  let high = Math.max(workload.chance, candidate.asimBaseline - minimumLoss);
  if (high === workload.chance) high = Math.min(candidate.asimBaseline, workload.chance + (workload.task === "ImageNet" ? 1 : 5));
  const insideAdc = adcDeficit === 0;
  const insideNoise = noiseRatio <= 1;
  const evidenceClass = insideAdc && insideNoise ? "Published low-loss regime" : (adcDeficit <= 1 && noiseRatio <= 2 ? "Boundary regime" : "Outside low-loss regime");
  const encodingText = encoding === 0 ? "None" : String(encoding);
  const configLines = [
    `# Screening translation; paste into the matching ASiM example config`,
    `cfg.task = '${workload.task}'`,
    `# model = ${workload.model}; checkpoint baseline = ${candidate.asimBaseline.toFixed(1)}% top-1`
  ];
  if (workload.model === "ResNet-18") {
    ["conv", "linear"].forEach((module) => {
      configLines.push(
        `cfg.asim_${module}_wbit = ${candidate.weightBits}`,
        `cfg.asim_${module}_xbit = ${candidate.inputBits}`,
        `cfg.asim_${module}_adc_prec = ${candidate.adcBits}`,
        `cfg.asim_${module}_nrow = ${Math.round(candidate.dimension)}`,
        `cfg.asim_${module}_rand_noise_sigma = ${candidate.asimRandomNoise}`,
        `cfg.asim_${module}_non_linear_sigma = ${candidate.asimNonlinearity}`,
        `cfg.asim_${module === "conv" ? "act_enc" : "linear_act_enc"} = ${encodingText}`
      );
    });
  } else {
    configLines.push(
      `cfg.asim_vit_attn_nrow = ${Math.round(candidate.dimension)}`,
      `cfg.asim_vit_attn_qk_qbit = ${candidate.inputBits}`,
      `cfg.asim_vit_attn_qk_kbit = ${candidate.inputBits}`,
      `cfg.asim_vit_attn_av_abit = ${candidate.inputBits}`,
      `cfg.asim_vit_attn_av_vbit = ${candidate.inputBits}`,
      `cfg.asim_vit_attn_qk_adc_prec = ${candidate.adcBits}`,
      `cfg.asim_vit_attn_av_adc_prec = ${candidate.adcBits}`,
      `cfg.asim_vit_attn_proj_adc_prec = ${candidate.adcBits}`,
      `cfg.asim_vit_attn_proj_wbit = ${candidate.weightBits}`,
      `cfg.asim_vit_attn_proj_xbit = ${candidate.inputBits}`,
      `cfg.asim_vit_attn_qk_rand_noise_sigma = ${candidate.asimRandomNoise}`,
      `cfg.asim_vit_attn_av_rand_noise_sigma = ${candidate.asimRandomNoise}`,
      `cfg.asim_vit_attn_proj_rand_noise_sigma = ${candidate.asimRandomNoise}`,
      `cfg.asim_vit_attn_qk_non_linear_sigma = ${candidate.asimNonlinearity}`,
      `cfg.asim_vit_attn_av_non_linear_sigma = ${candidate.asimNonlinearity}`,
      `cfg.asim_vit_attn_proj_non_linear_sigma = ${candidate.asimNonlinearity}`,
      `cfg.asim_vit_attn_qk_k_enc = ${encodingText}`,
      `cfg.asim_vit_attn_av_a_enc = ${encodingText}`,
      `cfg.asim_vit_attn_proj_act_enc = ${encodingText}`
    );
    ["mlp", "fc"].forEach((module) => {
      configLines.push(
        `cfg.asim_vit_${module}_wbit = ${candidate.weightBits}`,
        `cfg.asim_vit_${module}_xbit = ${candidate.inputBits}`,
        `cfg.asim_vit_${module}_adc_prec = ${candidate.adcBits}`,
        `cfg.asim_vit_${module}_nrow = ${Math.round(candidate.dimension)}`,
        `cfg.asim_vit_${module}_rand_noise_sigma = ${candidate.asimRandomNoise}`,
        `cfg.asim_vit_${module}_non_linear_sigma = ${candidate.asimNonlinearity}`,
        `cfg.asim_vit_${module}_act_enc = ${encodingText}`
      );
    });
    configLines.push(
      `cfg.asim_vit_quant_conv_wbit = ${candidate.weightBits}`,
      `cfg.asim_vit_quant_conv_xbit = ${candidate.inputBits}`
    );
  }
  configLines.push(
    `# Model preparation selected here: ${candidate.asimTraining === "nat" ? "noise-aware training (NAT)" : "quantization-aware training (QAT)"}.`,
    `# Run main/src_simulation.py to obtain the actual validation-set Test Acc.`
  );
  return { workload, baseAdc, fullRangeAdc, requiredAdc, adcDeficit, combinedNoisePct, noiseLsb, noiseRatio, low, high, evidenceClass, config: configLines.join("\n") };
}

async function initDetailedLab() {
  const chart = document.querySelector("#lab-chart");
  if (!chart) return;
  const rawRows = await loadBenchmarkRows();
  const rows = rawRows.map((row) => ({
    ...row,
    Architecture: row.Architecture.trim(),
    model: row["Compute Model"].trim(),
    tech: labNumber(row["Tech (nm)"]),
    dimension: labNumber(row.N),
    inputBits: labNumber(row.B_x),
    weightBits: labNumber(row.B_w),
    adcBits: labNumber(row.B_ADC),
    efficiency: labNumber(row["TOPS/W"]),
    density: labNumber(row["TOPS/mm2"]),
    throughput: labNumber(row.TOPS)
  }));
  const plottedRows = rows.filter((row) => row.efficiency > 0 && row.density > 0);

  const controls = {
    architecture: document.querySelector("#lab-architecture"), model: document.querySelector("#lab-model"),
    tech: document.querySelector("#lab-tech"), dimension: document.querySelector("#lab-dimension"),
    inputBits: document.querySelector("#lab-input-bits"), weightBits: document.querySelector("#lab-weight-bits"),
    energy: document.querySelector("#lab-energy"), throughput: document.querySelector("#lab-throughput"), area: document.querySelector("#lab-area"),
    adcBits: document.querySelector("#lab-adc-bits"), information: document.querySelector("#lab-information"),
    hasSnr: document.querySelector("#lab-has-snr"), snr: document.querySelector("#lab-snr"),
    envmDevice: document.querySelector("#lab-envm-device"), vbl: document.querySelector("#lab-vbl"), corePeriod: document.querySelector("#lab-core-period"),
    cellCap: document.querySelector("#lab-cell-cap"), adcNoise: document.querySelector("#lab-adc-noise"),
    asimWorkload: document.querySelector("#lab-asim-workload"), asimBaseline: document.querySelector("#lab-asim-baseline"),
    asimEncoding: document.querySelector("#lab-asim-encoding"), asimRandomNoise: document.querySelector("#lab-asim-random-noise"),
    asimNonlinearity: document.querySelector("#lab-asim-nonlinearity"), asimTraining: document.querySelector("#lab-asim-training")
  };
  const outputs = {
    tech: document.querySelector("#lab-tech-output"), energy: document.querySelector("#lab-energy-output"),
    throughput: document.querySelector("#lab-throughput-output"), area: document.querySelector("#lab-area-output"),
    information: document.querySelector("#lab-information-output"), snr: document.querySelector("#lab-snr-output"),
    vbl: document.querySelector("#lab-vbl-output"), corePeriod: document.querySelector("#lab-core-period-output"),
    cellCap: document.querySelector("#lab-cell-cap-output"), adcNoise: document.querySelector("#lab-adc-noise-output"),
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

  function candidateValues() {
    const energyFj = 10 ** Number(controls.energy.value);
    const throughput = 10 ** Number(controls.throughput.value);
    const area = 10 ** Number(controls.area.value);
    return {
      architecture: controls.architecture.value, model: controls.model.value,
      tech: Number(controls.tech.value), dimension: Number(controls.dimension.value),
      inputBits: Number(controls.inputBits.value), weightBits: Number(controls.weightBits.value),
      adcBits: controls.adcBits.value ? Number(controls.adcBits.value) : null,
      information: Number(controls.information.value), energyFj, throughput, area,
      efficiency: 1000 / energyFj, density: throughput / area,
      hasSnr: controls.hasSnr.checked, snr: Number(controls.snr.value),
      envmDevice: controls.envmDevice.value, vbl: 10 ** Number(controls.vbl.value), corePeriodNs: Number(controls.corePeriod.value),
      cellCap: Number(controls.cellCap.value), adcNoiseMv: Number(controls.adcNoise.value),
      asimWorkload: controls.asimWorkload.value, asimBaseline: Number(controls.asimBaseline.value),
      asimEncoding: Number(controls.asimEncoding.value), asimRandomNoise: Number(controls.asimRandomNoise.value),
      asimNonlinearity: Number(controls.asimNonlinearity.value), asimTraining: controls.asimTraining.value
    };
  }

  function refreshModelOptions() {
    const architecture = controls.architecture.value;
    const supported = new Set(rows.filter((row) => row.Architecture === architecture && row.model).map((row) => row.model));
    if (architecture === "eNVM") {
      supported.clear();
      supported.add("IS");
    }
    const prior = controls.model.value;
    controls.model.replaceChildren(...labModelOrder.filter((model) => supported.has(model)).map((model) => new Option(labModelLabels[model] || model, model)));
    if ([...supported].includes(prior)) controls.model.value = prior;
    refreshModelScope();
  }

  function refreshModelScope() {
    const architecture = controls.architecture.value;
    const model = controls.model.value;
    const modelRows = rows.filter((row) => row.Architecture === architecture && row.model === model);
    const dimensions = [...new Set(modelRows.map((row) => row.dimension).filter(Number.isFinite))].sort((a, b) => a - b);
    document.querySelector("#lab-dimension-help").textContent = dimensions.length ? `Reported for this pair: ${dimensions.join(", ")}.` : "No reported dimension for this pair.";
    const digitalImc = model === "DIMC";
    controls.adcBits.disabled = digitalImc;
    if (digitalImc) controls.adcBits.value = "";
    else if (!controls.adcBits.value) controls.adcBits.value = modelRows.some((row) => row.adcBits === 6) ? "6" : "5";
    document.querySelector("#envm-model-controls").hidden = !(architecture === "eNVM" && model === "IS");
    document.querySelector("#sram-model-controls").hidden = !(architecture === "SRAM" && model === "QR");
    document.querySelector("#asim-model-controls").hidden = !(architecture === "SRAM" && ["QS", "QR", "QS-QR"].includes(model));
  }

  function updateEvidence(candidate) {
    const pair = rows.filter((row) => row.Architecture === candidate.architecture && row.model === candidate.model);
    const node = pair.filter((row) => labEqual(row.tech, candidate.tech));
    const dimension = node.filter((row) => labEqual(row.dimension, candidate.dimension));
    const precision = dimension.filter((row) => labEqual(row.inputBits, candidate.inputBits) && labEqual(row.weightBits, candidate.weightBits));
    const exact = precision.filter((row) => candidate.adcBits === null ? row.adcBits === null : labEqual(row.adcBits, candidate.adcBits));
    const sameExceptNode = pair.filter((row) => labEqual(row.dimension, candidate.dimension) && labEqual(row.inputBits, candidate.inputBits) && labEqual(row.weightBits, candidate.weightBits) && (candidate.adcBits === null ? row.adcBits === null : labEqual(row.adcBits, candidate.adcBits)));
    const status = document.querySelector("#evidence-status");
    const trail = `${pair.length} model point${pair.length === 1 ? "" : "s"} → ${node.length} at ${candidate.tech} nm → ${dimension.length} at N=${candidate.dimension} → ${precision.length} at ${candidate.inputBits}b×${candidate.weightBits}b → ${exact.length} exact`;
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

  function validateModel(candidate) {
    if (candidate.architecture === "eNVM" && candidate.model === "IS" && window.IMC_MODEL_DATA) {
      const source = window.IMC_MODEL_DATA;
      const device = source.devices[candidate.envmDevice];
      const rowBracket = labBracket(source.dimensions, candidate.dimension, true);
      const columnBracket = labBracket(device.vbl, candidate.vbl, true);
      const modelDimension = Math.max(source.dimensions[0], Math.min(source.dimensions[source.dimensions.length - 1], candidate.dimension));
      const geq = labBilinear(device.averageGeq, rowBracket, columnBracket);
      const sndrDb = labBilinear(device.sndrDb, rowBracket, columnBracket);
      const adcBits = candidate.adcBits || 6;
      const periodSeconds = candidate.corePeriodNs * 1e-9;
      const constants = source.energyConstants;
      const arrayJ = device.vdd * candidate.vbl * periodSeconds * geq / (2 * modelDimension);
      const adcJ = (constants.adcLinearJPerBit * adcBits + constants.adcExponentialJ * 4 ** adcBits) / (2 * modelDimension);
      const energyFj = (arrayJ + adcJ) * 1e15;
      const perCoreTops = (2 * modelDimension / periodSeconds) / 1e12;
      const withinSurface = !rowBracket.clamped && !columnBracket.clamped;
      return {
        lead: candidate.energyFj + 1e-9 >= energyFj ? "Claim is not below the modeled core-energy floor" : "Claim is below the modeled core-energy floor",
        detail: `${candidate.envmDevice} current-summing model at N=${candidate.dimension}, ${labFormat(candidate.vbl, 4)} V, and ${candidate.corePeriodNs} ns. The source surface uses a six-bit ADC and omits system-level overhead${withinSurface ? "." : "; the selected point was clamped to its published parameter range."}`,
        metrics: [
          ["Modeled core energy", `${labFormat(energyFj, 3)} fJ/1b-op`],
          ["Modeled efficiency ceiling", `${labFormat(1000 / energyFj, 1)} 1b-TOPS/W`],
          ["Per-core throughput", `${labFormat(perCoreTops, 4)} 1b-TOPS`],
          ["Implied parallel cores", labFormat(candidate.throughput / perCoreTops, 1)],
          ["Published-surface SNDR", `${labFormat(sndrDb, 1)} dB${candidate.adcBits === 6 ? "" : " (6-bit reference)"}`]
        ],
        modeledSnr: candidate.adcBits === 6 ? sndrDb : null,
        referenceSnr: sndrDb,
        referenceNode: "22 nm MRAM-validated behavioral parameters; no automatic node scaling"
      };
    }
    if (candidate.architecture === "SRAM" && candidate.model === "QR") {
      const sndrDb = labSramQrSnr(candidate);
      const atSourceNode = candidate.tech === 28;
      return {
        lead: sndrDb === null ? "Selected precision is outside the QR simulator" : atSourceNode ? "28 nm QR compute-SNR estimate available" : "Only a 28 nm QR reference is available",
        detail: sndrDb === null ? "The SRAM repository simulator requires integer input, weight, and ADC precision." : `Deterministic capacitance-mismatch and ADC-noise simulation gives ${labFormat(sndrDb, 1)} dB at the repository's 28 nm parameterization. The repository does not model core energy, area, or throughput.`,
        metrics: sndrDb === null ? [] : [
          ["28 nm reference SNR", `${labFormat(sndrDb, 1)} dB`],
          ["Unit capacitance", `${labFormat(candidate.cellCap, 2)} fF`],
          ["ADC noise", `${labFormat(candidate.adcNoiseMv, 2)} mV`],
          ["Energy / area validation", "Not provided by source model"]
        ],
        modeledSnr: atSourceNode ? sndrDb : null,
        referenceSnr: sndrDb,
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
    title.textContent = asimActive ? "Network-level accuracy" : "Accuracy awareness";
    scale.hidden = asimActive;
    asimResult.hidden = !asimActive;
    if (asimActive) {
      marker.hidden = true;
      if (!Number.isInteger(candidate.inputBits) || !Number.isInteger(candidate.weightBits) || candidate.adcBits === null) {
        lead.textContent = "ASiM screening needs integer precisions and an ADC";
        detail.textContent = "Choose integer input and weight bits plus an ADC precision. ASiM decomposes quantized tensors into binary MAC cycles and cannot use the selected fractional precision directly.";
        document.querySelector("#asim-adc-boundary").textContent = "Not evaluated";
        document.querySelector("#asim-noise").textContent = "Not evaluated";
        document.querySelector("#asim-evidence-class").textContent = "Outside ASiM input scope";
        document.querySelector("#asim-config-text").textContent = "Select integer input/weight precision and an ADC precision to generate an ASiM configuration.";
        return;
      }
      const screening = labAsimScreen(candidate);
      lead.textContent = `Screening band: ${labFormat(screening.low, 1)}–${labFormat(screening.high, 1)}% top-1`;
      const trainingNote = candidate.asimTraining === "nat" ? "The entered baseline should be the clean accuracy of the selected NAT checkpoint; its added robustness is not credited numerically without a direct run." : "The entered baseline should be the clean accuracy of the selected QAT checkpoint.";
      detail.textContent = `${screening.workload.label} screening band from ASiM’s published ADC and noise regimes, anchored to the entered ${labFormat(candidate.asimBaseline, 1)}% digital baseline. ${trainingNote} This is not an ASiM inference result.`;
      document.querySelector("#asim-adc-boundary").textContent = `${screening.requiredAdc}b guidance · ${screening.fullRangeAdc}b full range`;
      document.querySelector("#asim-noise").textContent = `${labFormat(screening.combinedNoisePct, 3)}% Vpp · ${labFormat(screening.noiseLsb, 3)} LSB rms proxy`;
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
    const pointGroup = svgElement("g", { "aria-hidden": "true" });
    plottedRows.forEach((row) => pointGroup.append(svgElement("circle", { cx: xPosition(row.density), cy: yPosition(row.efficiency), r: compact ? 5.2 : 4.2, fill: colors[row.Architecture] || "#687080", class: "reported-point" })));
    chart.append(pointGroup);
    const x = xPosition(Math.max(xDomain[0], Math.min(xDomain[1], candidate.density)));
    const y = yPosition(Math.max(yDomain[0], Math.min(yDomain[1], candidate.efficiency)));
    const size = compact ? 9 : 10;
    chart.append(svgElement("polygon", { points: `${x},${y-size} ${x+size},${y} ${x},${y+size} ${x-size},${y}`, fill: colors[candidate.architecture] || "#20242b", class: "candidate-point" }));
    const label = svgElement("text", { x: Math.min(x + 15, width - 118), y: Math.max(y - 14, 24), class: "candidate-label" }); label.textContent = `Your ${candidate.architecture}`; chart.append(label);
  }

  function updateNeighbors(candidate) {
    const nearestByIndex = new Map();
    plottedRows.forEach((row) => {
      const distance = Math.hypot(Math.log10(row.density / candidate.density), Math.log10(row.efficiency / candidate.efficiency));
      const prior = nearestByIndex.get(row.Index);
      if (!prior || distance < prior.distance) nearestByIndex.set(row.Index, { ...row, distance });
    });
    const nearest = [...nearestByIndex.values()].sort((a, b) => a.distance - b.distance).slice(0, 3);
    document.querySelector("#neighbor-body").innerHTML = nearest.map((row) => `<tr><td><span class="index-chip">#${escapeHtml(row.Index)}</span></td><td>${escapeHtml(row["Paper Title"])}</td><td><span class="arch-chip ${escapeHtml(row.Architecture.toLowerCase())}">${escapeHtml(row.Architecture)}</span></td></tr>`).join("");
  }

  function update() {
    const candidate = candidateValues();
    outputs.tech.textContent = `${candidate.tech} nm`;
    outputs.energy.textContent = `${labFormat(candidate.energyFj, 3)} fJ`;
    outputs.throughput.textContent = `${labFormat(candidate.throughput, 3)} TOPS`;
    outputs.area.textContent = `${labFormat(candidate.area, 4)} mm²`;
    outputs.information.textContent = `${labFormat(candidate.information, 1)} bits`;
    outputs.snr.textContent = `${labFormat(candidate.snr, 1)} dB`;
    outputs.vbl.textContent = `${labFormat(candidate.vbl, 4)} V`;
    outputs.corePeriod.textContent = `${candidate.corePeriodNs} ns`;
    outputs.cellCap.textContent = `${labFormat(candidate.cellCap, 2)} fF`;
    outputs.adcNoise.textContent = `${labFormat(candidate.adcNoiseMv, 2)} mV`;
    outputs.asimRandomNoise.innerHTML = `${labFormat(candidate.asimRandomNoise, 3)}% V<sub>pp</sub>`;
    outputs.asimNonlinearity.innerHTML = `${labFormat(candidate.asimNonlinearity, 3)}% V<sub>pp</sub>`;
    controls.snr.disabled = !candidate.hasSnr;
    if (candidate.envmDevice === "FeFET") controls.vbl.max = "0.7"; else controls.vbl.max = "-0.09";

    const evidence = updateEvidence(candidate);
    const modelResult = validateModel(candidate);
    draw(candidate);
    document.querySelector("#placement-lead").textContent = `${candidate.architecture} · ${candidate.model} · ${candidate.tech} nm`;
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
    pendingUpdate = setTimeout(update, 70);
  }
  controls.architecture.addEventListener("change", () => { refreshModelOptions(); update(); });
  controls.model.addEventListener("change", () => { refreshModelScope(); update(); });
  controls.envmDevice.addEventListener("change", () => {
    const maximum = controls.envmDevice.value === "FeFET" ? 0.7 : -0.09;
    controls.vbl.max = String(maximum);
    if (Number(controls.vbl.value) > maximum) controls.vbl.value = String(maximum);
    update();
  });
  const asimBaselineDefaults = { "resnet18-cifar10": 93, "resnet18-imagenet": 69.8, "vitb32-cifar10": 98, "vitb32-imagenet": 75.9 };
  controls.asimWorkload.addEventListener("change", () => {
    controls.asimBaseline.value = String(asimBaselineDefaults[controls.asimWorkload.value]);
    update();
  });
  Object.entries(controls).filter(([key]) => !["architecture", "model", "envmDevice"].includes(key)).forEach(([, control]) => control.addEventListener("input", scheduleUpdate));
  refreshModelOptions();
  update();
}

initDetailedLab().catch((error) => {
  const lead = document.querySelector("#placement-lead");
  if (lead) lead.textContent = error.message;
});
