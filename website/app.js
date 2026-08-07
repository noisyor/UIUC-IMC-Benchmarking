const menuButton = document.querySelector(".menu-toggle");
const siteNav = document.querySelector("#site-nav");

if (menuButton && siteNav) {
  menuButton.addEventListener("click", () => {
    const open = menuButton.getAttribute("aria-expanded") !== "true";
    menuButton.setAttribute("aria-expanded", String(open));
    menuButton.textContent = open ? "Close" : "Menu";
    siteNav.classList.toggle("open", open);
    document.body.classList.toggle("menu-open", open);
  });
  siteNav.querySelectorAll("a").forEach((link) => link.addEventListener("click", () => {
    menuButton.setAttribute("aria-expanded", "false");
    menuButton.textContent = "Menu";
    siteNav.classList.remove("open");
    document.body.classList.remove("menu-open");
  }));
}

const plotLayouts = {
  efficiency: [
    [11,77],[23,61],[32,71],[40,48],[49,58],[57,35],[66,43],[75,24],[83,31],[88,14],[69,68],[44,80]
  ],
  throughput: [
    [8,84],[19,70],[27,55],[36,67],[46,40],[54,52],[64,27],[73,43],[80,20],[90,34],[59,73],[41,24]
  ],
  precision: [
    [10,30],[20,43],[29,25],[38,58],[47,38],[56,72],[64,48],[74,65],[83,42],[90,79],[69,22],[43,82]
  ]
};

document.querySelectorAll("[data-axis]").forEach((button) => {
  button.addEventListener("click", () => {
    const layout = plotLayouts[button.dataset.axis];
    document.querySelectorAll("[data-axis]").forEach((item) => item.setAttribute("aria-pressed", String(item === button)));
    document.querySelectorAll(".plot-point").forEach((point, index) => {
      point.style.left = `${layout[index][0]}%`;
      point.style.bottom = `${layout[index][1]}%`;
    });
    const labels = {
      efficiency: ["1b-TOPS/mm²", "1b-TOPS/W"],
      throughput: ["1b-TOPS", "1b-TOPS/W"],
      precision: ["Pre-ADC information", "ADC precision"]
    };
    document.querySelector(".plot-axis-x").textContent = labels[button.dataset.axis][0];
    document.querySelector(".plot-axis-y").textContent = labels[button.dataset.axis][1];
  });
});

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let index = 0; index < text.length; index += 1) {
    const char = text[index];
    if (char === '"') {
      if (quoted && text[index + 1] === '"') { field += '"'; index += 1; }
      else quoted = !quoted;
    } else if (char === "," && !quoted) {
      row.push(field); field = "";
    } else if ((char === "\n" || char === "\r") && !quoted) {
      if (char === "\r" && text[index + 1] === "\n") index += 1;
      row.push(field); field = "";
      if (row.some((value) => value.trim())) rows.push(row);
      row = [];
    } else field += char;
  }
  if (field || row.length) { row.push(field); rows.push(row); }
  const headers = rows.shift().map((header) => header.trim());
  return rows.map((values) => Object.fromEntries(headers.map((header, index) => [header, (values[index] || "").trim()])));
}

async function loadBenchmarkRows() {
  if (typeof window.BENCHMARK_CSV === "string") return parseCsv(window.BENCHMARK_CSV);
  const response = await fetch("data/Benchmarking_Data.csv");
  if (!response.ok) throw new Error("Benchmark data could not be loaded.");
  return parseCsv(await response.text());
}

function titleCaseArch(value) {
  const labels = { SRAM: "SRAM", eNVM: "eNVM", eDRAM: "eDRAM", eFlash: "eFlash", Digital: "Digital" };
  return labels[value] || value;
}

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (character) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;"
  })[character]);
}

function paperPeople(index, fallbackAuthors) {
  const metadata = window.BENCHMARK_PAPER_METADATA?.papers?.[index];
  const authors = metadata?.authors || [];
  if (!authors.length) return {
    toggle: `<span class="paper-meta">${escapeHtml(fallbackAuthors)}</span>`,
    row: ""
  };
  const affiliations = [...new Set(authors.flatMap((author) => author.affiliations || []))];
  const people = authors.map((author) => {
    const authorAffiliations = author.affiliations?.length
      ? author.affiliations.map((affiliation) => escapeHtml(affiliation)).join("<br>")
      : "Affiliation unavailable in the DOI metadata";
    return `<div class="paper-person"><dt>${escapeHtml(author.name)}</dt><dd>${authorAffiliations}</dd></div>`;
  }).join("");
  const sourceLink = metadata.openAlexUrl
    ? `<a href="${escapeHtml(metadata.openAlexUrl)}" target="_blank" rel="noopener">OpenAlex record ↗</a>`
    : "OpenAlex DOI metadata";
  const label = `${authors.length} author${authors.length === 1 ? "" : "s"} · ${affiliations.length} affiliation${affiliations.length === 1 ? "" : "s"}`;
  const panelId = `paper-people-${index}`;
  return {
    toggle: `<button class="paper-people-toggle" type="button" aria-expanded="false" aria-controls="${panelId}">${label}</button>`,
    row: `
      <tr id="${panelId}" class="paper-people-row" hidden>
        <td colspan="9">
          <section class="paper-people-panel" aria-label="Authors and affiliations for benchmark paper ${escapeHtml(index)}">
            <div class="paper-people-head"><strong>Full author and affiliation record</strong><span>Source: ${sourceLink}</span></div>
            <dl class="paper-people-list">${people}</dl>
          </section>
        </td>
      </tr>`
  };
}

async function initExplorer() {
  const tableBody = document.querySelector("#benchmark-body");
  if (!tableBody) return;
  const operatingPoints = await loadBenchmarkRows();
  const papers = new Map();
  operatingPoints.forEach((row) => {
    if (!row.Index) return;
    if (!papers.has(row.Index)) papers.set(row.Index, { ...row, operatingPoints: 0, ranges: {} });
    const paper = papers.get(row.Index);
    paper.operatingPoints += 1;
    ["TOPS/W", "TOPS/mm2", "TOPS"].forEach((metric) => {
      const value = Number(row[metric]);
      if (!Number.isFinite(value)) return;
      if (!paper.ranges[metric]) paper.ranges[metric] = [value, value];
      paper.ranges[metric][0] = Math.min(paper.ranges[metric][0], value);
      paper.ranges[metric][1] = Math.max(paper.ranges[metric][1], value);
    });
  });
  const paperRows = [...papers.values()].sort((a, b) => Number(b.Index) - Number(a.Index));
  const search = document.querySelector("#paper-search");
  const architecture = document.querySelector("#architecture-filter");
  const year = document.querySelector("#year-filter");
  const venue = document.querySelector("#venue-filter");
  const count = document.querySelector("#result-count");

  tableBody.addEventListener("click", (event) => {
    const toggle = event.target.closest(".paper-people-toggle");
    if (!toggle) return;
    const detailRow = document.querySelector(`#${toggle.getAttribute("aria-controls")}`);
    if (!detailRow) return;
    const expanded = toggle.getAttribute("aria-expanded") === "true";
    toggle.setAttribute("aria-expanded", String(!expanded));
    detailRow.hidden = expanded;
  });

  [...new Set(paperRows.map((row) => row.Year))].sort().reverse().forEach((value) => year.add(new Option(value, value)));
  [...new Set(paperRows.map((row) => row["Journal/Conference"]))].sort().forEach((value) => venue.add(new Option(value, value)));

  const metric = (paper, name) => {
    const range = paper.ranges[name];
    if (!range) return "—";
    const format = (value) => Number(value).toLocaleString(undefined, { maximumFractionDigits: 2 });
    return range[0] === range[1] ? format(range[0]) : `${format(range[0])}–${format(range[1])}`;
  };

  function render() {
    const query = search.value.trim().toLowerCase();
    const visible = paperRows.filter((row) => {
      const metadata = window.BENCHMARK_PAPER_METADATA?.papers?.[row.Index];
      const people = metadata?.authors?.flatMap((author) => [author.name, ...(author.affiliations || [])]) || [];
      const haystack = [row.Index, row["Paper Title"], row.Authors, row["Journal/Conference"], ...people].join(" ").toLowerCase();
      return (!query || haystack.includes(query)) &&
        (!architecture.value || row.Architecture === architecture.value) &&
        (!year.value || row.Year === year.value) &&
        (!venue.value || row["Journal/Conference"] === venue.value);
    });
    count.textContent = `${visible.length} indexed paper${visible.length === 1 ? "" : "s"}`;
    if (!visible.length) {
      tableBody.innerHTML = '<tr><td class="explorer-empty" colspan="9">No papers match these filters.</td></tr>';
      return;
    }
    tableBody.innerHTML = visible.map((row) => {
      const paperLink = window.BENCHMARK_PAPER_LINKS?.[row.Index];
      const people = paperPeople(row.Index, row.Authors);
      const linkedTitle = paperLink
        ? `<a class="paper-title-link" href="${escapeHtml(paperLink.url)}" target="_blank" rel="noopener">${escapeHtml(row["Paper Title"])}</a>`
        : escapeHtml(row["Paper Title"]);
      const linkCell = paperLink
        ? `<a class="paper-link" href="${escapeHtml(paperLink.url)}" target="_blank" rel="noopener" aria-label="Open DOI record for ${escapeHtml(row["Paper Title"])}">DOI ↗</a>`
        : "—";
      return `
        <tr>
          <td><span class="index-chip">#${escapeHtml(row.Index)}</span></td>
          <td>${escapeHtml(row.Year)}</td>
          <td><strong>${linkedTitle}</strong>${people.toggle}</td>
          <td>${escapeHtml(row["Journal/Conference"])}</td>
          <td><span class="arch-chip ${escapeHtml(row.Architecture.toLowerCase())}">${escapeHtml(titleCaseArch(row.Architecture))}</span></td>
          <td>${metric(row, "TOPS/W")}</td>
          <td>${metric(row, "TOPS/mm2")}</td>
          <td>${row.operatingPoints}</td>
          <td>${linkCell}</td>
        </tr>${people.row}`;
    }).join("");
  }

  [search, architecture, year, venue].forEach((control) => control.addEventListener(control === search ? "input" : "change", render));
  render();
}

initExplorer().catch((error) => {
  const tableBody = document.querySelector("#benchmark-body");
  if (tableBody) tableBody.innerHTML = `<tr><td class="explorer-empty" colspan="8">${error.message} Serve the site over HTTP to use the explorer.</td></tr>`;
});

function svgElement(name, attributes = {}) {
  const element = document.createElementNS("http://www.w3.org/2000/svg", name);
  Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
  return element;
}

function percentile(values, candidate) {
  if (!values.length) return 0;
  return Math.round(100 * values.filter((value) => value <= candidate).length / values.length);
}

async function initLab() {
  const chart = document.querySelector("#lab-chart");
  if (!chart) return;
  const allRows = (await loadBenchmarkRows()).map((row) => ({
    ...row,
    efficiency: Number(row["TOPS/W"]),
    density: Number(row["TOPS/mm2"])
  })).filter((row) => row.efficiency > 0 && row.density > 0);

  const controls = {
    architecture: document.querySelector("#lab-architecture"),
    model: document.querySelector("#lab-model"),
    tech: document.querySelector("#lab-tech"),
    dimension: document.querySelector("#lab-dimension"),
    inputBits: document.querySelector("#lab-input-bits"),
    weightBits: document.querySelector("#lab-weight-bits"),
    efficiency: document.querySelector("#lab-efficiency"),
    density: document.querySelector("#lab-density"),
    adcBits: document.querySelector("#lab-adc-bits"),
    information: document.querySelector("#lab-information"),
    hasSnr: document.querySelector("#lab-has-snr"),
    snr: document.querySelector("#lab-snr")
  };
  const outputs = {
    tech: document.querySelector("#lab-tech-output"),
    efficiency: document.querySelector("#lab-efficiency-output"),
    density: document.querySelector("#lab-density-output"),
    information: document.querySelector("#lab-information-output"),
    snr: document.querySelector("#lab-snr-output")
  };
  const colors = { SRAM: "#218a5b", eNVM: "#df4b45", eFlash: "#df4b45", Digital: "#3268cc", eDRAM: "#d98a1c" };
  const measuredWidth = Math.round(chart.getBoundingClientRect().width || 900);
  const width = Math.max(320, measuredWidth);
  const compactChart = width < 560;
  const height = compactChart ? 420 : 560;
  const margin = compactChart ? { top: 20, right: 16, bottom: 66, left: 72 } : { top: 24, right: 28, bottom: 78, left: 96 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const xDomain = [0.003, 6000];
  const yDomain = [4, 25000];
  const xTicks = compactChart ? [0.01, 1, 100, 1000] : [0.01, 0.1, 1, 10, 100, 1000];
  const yTicks = [10, 100, 1000, 10000];
  const logPosition = (value, domain, size) => (Math.log10(value) - Math.log10(domain[0])) / (Math.log10(domain[1]) - Math.log10(domain[0])) * size;
  const xPosition = (value) => margin.left + logPosition(value, xDomain, plotWidth);
  const yPosition = (value) => margin.top + plotHeight - logPosition(value, yDomain, plotHeight);
  const formatMetric = (value) => value.toLocaleString(undefined, { maximumFractionDigits: value < 1 ? 3 : 1 });
  chart.setAttribute("viewBox", `0 0 ${width} ${height}`);

  function drawBase() {
    chart.replaceChildren(chart.querySelector("title"), chart.querySelector("desc"));
    const grid = svgElement("g");
    xTicks.forEach((tick) => {
      const x = xPosition(tick);
      grid.append(svgElement("line", { x1: x, y1: margin.top, x2: x, y2: margin.top + plotHeight, class: "grid-line" }));
      const label = svgElement("text", { x, y: margin.top + plotHeight + 28, "text-anchor": "middle" });
      label.textContent = tick;
      grid.append(label);
    });
    yTicks.forEach((tick) => {
      const y = yPosition(tick);
      grid.append(svgElement("line", { x1: margin.left, y1: y, x2: margin.left + plotWidth, y2: y, class: "grid-line" }));
      const label = svgElement("text", { x: margin.left - 14, y: y + 5, "text-anchor": "end" });
      label.textContent = tick.toLocaleString();
      grid.append(label);
    });
    grid.append(svgElement("line", { x1: margin.left, y1: margin.top + plotHeight, x2: margin.left + plotWidth, y2: margin.top + plotHeight, class: "axis-line" }));
    grid.append(svgElement("line", { x1: margin.left, y1: margin.top, x2: margin.left, y2: margin.top + plotHeight, class: "axis-line" }));
    const xTitle = svgElement("text", { x: margin.left + plotWidth / 2, y: height - 18, "text-anchor": "middle", class: "axis-title" });
    xTitle.textContent = "Compute density (1b-TOPS/mm²)";
    grid.append(xTitle);
    const yTitle = svgElement("text", { x: 24, y: margin.top + plotHeight / 2, "text-anchor": "middle", class: "axis-title", transform: `rotate(-90 24 ${margin.top + plotHeight / 2})` });
    yTitle.textContent = "Energy efficiency (1b-TOPS/W)";
    grid.append(yTitle);
    chart.append(grid);

    const points = svgElement("g", { "aria-hidden": "true" });
    allRows.forEach((row) => points.append(svgElement("circle", {
      cx: xPosition(row.density), cy: yPosition(row.efficiency), r: compactChart ? 5.2 : 4.2,
      fill: colors[row.Architecture] || "#687080", class: "reported-point"
    })));
    chart.append(points);
  }

  function candidateValues() {
    return {
      architecture: controls.architecture.value,
      model: controls.model.value,
      tech: Number(controls.tech.value),
      dimension: Number(controls.dimension.value),
      inputBits: Number(controls.inputBits.value),
      weightBits: Number(controls.weightBits.value),
      efficiency: 10 ** Number(controls.efficiency.value),
      density: 10 ** Number(controls.density.value),
      adcBits: Number(controls.adcBits.value),
      information: Number(controls.information.value),
      hasSnr: controls.hasSnr.checked,
      snr: Number(controls.snr.value)
    };
  }

  function updateAccuracy(candidate) {
    const lead = document.querySelector("#accuracy-lead");
    const detail = document.querySelector("#accuracy-detail");
    const marker = document.querySelector("#accuracy-marker");
    const gap = candidate.adcBits - candidate.information;
    if (!candidate.hasSnr) {
      lead.textContent = "Accuracy not established";
      detail.textContent = `ADC precision is ${Math.abs(gap).toFixed(1)} bits ${gap < 0 ? "below" : "above"} the selected pre-ADC information proxy. This comparison flags quantization pressure but cannot replace a measured compute SNR.`;
      marker.hidden = true;
      return;
    }
    marker.hidden = false;
    marker.style.left = `${Math.max(0, Math.min(100, candidate.snr / 50 * 100))}%`;
    if (candidate.snr < 10) {
      lead.textContent = "Below the cited guidance range";
      detail.textContent = `${candidate.snr.toFixed(1)} dB is below the 10–40 dB task-dependent range cited in the OJ-SSCS paper. Network accuracy still requires direct evaluation.`;
    } else if (candidate.snr <= 40) {
      lead.textContent = "Within the cited guidance range";
      detail.textContent = `${candidate.snr.toFixed(1)} dB falls inside the paper’s task-dependent 10–40 dB range. This is context, not a network-accuracy guarantee.`;
    } else {
      lead.textContent = "Above the cited guidance range";
      detail.textContent = `${candidate.snr.toFixed(1)} dB exceeds the cited range, but analog distortion, mapping, and the target network still need direct accuracy measurement.`;
    }
  }

  function updateNeighbors(candidate) {
    const nearestByIndex = new Map();
    allRows.forEach((row) => {
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
    outputs.efficiency.textContent = `${formatMetric(candidate.efficiency)} 1b-TOPS/W`;
    outputs.density.textContent = `${formatMetric(candidate.density)} 1b-TOPS/mm²`;
    outputs.information.textContent = `${candidate.information.toFixed(candidate.information % 1 ? 1 : 0)} bits`;
    outputs.snr.textContent = `${candidate.snr.toFixed(candidate.snr % 1 ? 1 : 0)} dB`;
    controls.snr.disabled = !candidate.hasSnr;

    drawBase();
    const x = xPosition(Math.max(xDomain[0], Math.min(xDomain[1], candidate.density)));
    const y = yPosition(Math.max(yDomain[0], Math.min(yDomain[1], candidate.efficiency)));
    const size = compactChart ? 9 : 10;
    chart.append(svgElement("polygon", {
      points: `${x},${y-size} ${x+size},${y} ${x},${y+size} ${x-size},${y}`,
      fill: colors[candidate.architecture] || "#20242b", class: "candidate-point"
    }));
    const label = svgElement("text", { x: Math.min(x + 15, width - 118), y: Math.max(y - 14, 24), class: "candidate-label" });
    label.textContent = `Your ${candidate.architecture}`;
    chart.append(label);

    const sameArchitecture = allRows.filter((row) => row.Architecture === candidate.architecture);
    const comparisonRows = sameArchitecture.length >= 3 ? sameArchitecture : allRows;
    document.querySelector("#placement-lead").textContent = `${candidate.architecture} · ${candidate.model} · ${candidate.tech} nm`;
    document.querySelector("#energy-per-op").textContent = `${(1000 / candidate.efficiency).toLocaleString(undefined, { maximumFractionDigits: 2 })} fJ`;
    document.querySelector("#efficiency-percentile").textContent = `${percentile(comparisonRows.map((row) => row.efficiency), candidate.efficiency)}th`;
    document.querySelector("#density-percentile").textContent = `${percentile(comparisonRows.map((row) => row.density), candidate.density)}th`;
    updateAccuracy(candidate);
    updateNeighbors(candidate);
  }

  Object.values(controls).forEach((control) => control.addEventListener("input", update));
  drawBase();
  update();
}
