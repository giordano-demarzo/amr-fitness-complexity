/* The Resistome Atlas — static dashboard front-end.
   Zero dependencies. Loads precomputed JSON from ./data and renders an
   interactive antibiotic co-resistance network, a per-pathogen resistance
   timeline, and a forward resistance-risk forecast. */
'use strict';

const SVGNS = 'http://www.w3.org/2000/svg';
const $ = s => document.querySelector(s);
const el = (tag, attrs = {}, kids = []) => {
  const ns = /^(svg|g|circle|line|text|rect|path|defs|filter|feGaussianBlur|feMerge|feMergeNode)$/.test(tag);
  const n = ns ? document.createElementNS(SVGNS, tag) : document.createElement(tag);
  for (const k in attrs) {
    if (k === 'class') n.setAttribute('class', attrs[k]);
    else if (k === 'text') n.textContent = attrs[k];
    else if (k === 'html') n.innerHTML = attrs[k];
    else n.setAttribute(k, attrs[k]);
  }
  (Array.isArray(kids) ? kids : [kids]).forEach(c => c && n.appendChild(c));
  return n;
};

const TIER = { 1: { c: '#2a9d8f', n: 'Access' }, 2: { c: '#e9c46a', n: 'Watch' }, 3: { c: '#c9184a', n: 'Reserve' }, 0: { c: '#5b6b7d', n: 'Unclassified' } };
const tierC = t => (TIER[t] || TIER[0]).c;

// magma-like scale for resistance fraction 0..1
const MAGMA = [[0, '#140b23'], [.2, '#3b0f70'], [.4, '#8c2981'], [.6, '#de4968'], [.8, '#fe9f6d'], [1, '#fcfdbf']];
function magma(t) {
  t = Math.max(0, Math.min(1, t));
  for (let i = 1; i < MAGMA.length; i++) {
    if (t <= MAGMA[i][0]) {
      const [a, ca] = MAGMA[i - 1], [b, cb] = MAGMA[i];
      return mix(ca, cb, (t - a) / (b - a));
    }
  }
  return MAGMA[MAGMA.length - 1][1];
}
function mix(h1, h2, u) {
  const p = h => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
  const a = p(h1), b = p(h2);
  const r = a.map((v, i) => Math.round(v + (b[i] - v) * u));
  return `rgb(${r[0]},${r[1]},${r[2]})`;
}
const italic = n => { const p = n.split(' '); return p.length > 1 ? `<i>${p[0]} ${p.slice(1).join(' ')}</i>` : `<i>${n}</i>`; };

const state = { meta: null, abx: [], abById: {}, net: null, paths: [], detail: null, sel: null, year: 2024, playing: null };

// ============================================================ boot
Promise.all(['meta', 'antibiotics', 'network', 'pathogens', 'pathogen_detail']
  .map(f => fetch(`data/${f}.json`).then(r => r.json())))
  .then(([meta, abx, net, paths, detail]) => {
    state.meta = meta; state.abx = abx; state.net = net; state.paths = paths; state.detail = detail;
    abx.forEach((a, i) => state.abById[i] = a);
    state.byName = {}; abx.forEach((a, i) => state.byName[a.name] = i);
    buildAbxIndex();
    initChrome(); buildNetwork(); initSearch(); initTime(); initModals();
    // shareable deep link: ?p=<pathogen>&y=<year>
    const q = new URLSearchParams(location.search);
    const pin = q.get('p');
    if (pin) {
      const hit = state.paths.find(x => x.name.toLowerCase() === pin.toLowerCase())
        || state.paths.find(x => x.name.toLowerCase().includes(pin.toLowerCase()));
      if (hit) {
        selectPathogen(hit.name);
        const yy = parseInt(q.get('y'), 10);
        if (yy && state.meta.years.includes(yy)) setYear(yy);
      }
    }
  })
  .catch(e => { document.body.innerHTML = `<p style="padding:40px;color:#f88">Failed to load data: ${e}. Serve this folder over HTTP (e.g. <code>python3 -m http.server</code>).</p>`; });

// ============================================================ chrome (kpis, legend, meta)
function initChrome() {
  const m = state.meta;
  const fmt = n => n.toLocaleString('en-US');
  const kpis = [
    { num: fmt(m.n_species), lab: 'species' },
    { num: fmt(m.n_antibiotics), lab: 'antibiotics' },
    { num: `${m.year_min}–${m.year_max}`, lab: 'surveillance' },
    { num: m.n_countries ? fmt(m.n_countries) : '—', lab: 'countries' },
    { num: fmt(m.n_records), lab: 'measurements' },
  ];
  $('#kpis').innerHTML = kpis.map(k =>
    `<div class="kpi"><div class="num">${k.num}</div><div class="lab">${k.lab}</div></div>`).join('');
  // method numbers
  const hubMax = Math.max(...state.abx.map(a => a.auc || 0));
  $('#m-auc').textContent = hubMax.toFixed(2);
  $('#m-rho').textContent = m.rho_complexity_inblock;
  // network legend
  $('#net-legend').innerHTML = [1, 2, 3].map(t =>
    `<span class="lg"><span class="dot" style="background:${tierC(t)}"></span>${TIER[t].n}</span>`).join('');
  // empty pathogen card (quick picks) — delegated so it survives re-renders
  renderEmptyCard();
  $('#pathcard').addEventListener('click', e => {
    const q = e.target.closest('.qp'); if (q) { selectPathogen(q.dataset.p); return; }
    if (e.target.closest('#openHeat')) $('#heatModal').hidden = false;
  });
  $('#deselect').addEventListener('click', deselectPathogen);
}

// ============================================================ modals (intro / about / timeline / antibiotic)
const closeAllModals = () => document.querySelectorAll('.modal-overlay').forEach(ov => ov.hidden = true);
function initModals() {
  const intro = $('#introModal'), about = $('#aboutModal');
  $('#openAbout').addEventListener('click', () => { about.hidden = false; });
  $('#introAbout').addEventListener('click', () => { intro.hidden = true; about.hidden = false; });
  // every [data-close] button closes; clicking a modal backdrop closes it; Esc closes all
  document.querySelectorAll('[data-close]').forEach(b => b.addEventListener('click', closeAllModals));
  document.querySelectorAll('.modal-overlay').forEach(ov =>
    ov.addEventListener('click', e => { if (e.target === ov) ov.hidden = true; }));
  document.addEventListener('keydown', e => { if (e.key === 'Escape') closeAllModals(); });
  // show the intro on first visit only (and not when arriving via a shared pathogen link)
  const pinned = new URLSearchParams(location.search).get('p');
  let seen = null; try { seen = localStorage.getItem('atlas_intro_seen'); } catch (e) {}
  if (!seen && !pinned) { intro.hidden = false; }
  try { localStorage.setItem('atlas_intro_seen', '1'); } catch (e) {}
}

// reverse indexes: for each antibiotic, which species already resist it, and which are most at risk next
function buildAbxIndex() {
  state.abxResist = {}; state.abxRisk = {};
  for (const pname in state.detail) {
    const d = state.detail[pname];
    for (const idx in d.lit) {
      const a = state.abById[idx]; if (!a) continue;
      const ts = d.timeline[a.name];
      const propR = ts ? ts[ts.length - 1][1] : null;
      (state.abxResist[a.name] = state.abxResist[a.name] || []).push({ p: pname, propR, year: d.lit[idx] });
    }
    (d.forecast || []).forEach(f => {
      (state.abxRisk[f.abx] = state.abxRisk[f.abx] || []).push({ p: pname, prob: f.prob, cur: f.cur_propR });
    });
  }
  for (const k in state.abxResist) state.abxResist[k].sort((a, b) => (b.propR || 0) - (a.propR || 0));
  for (const k in state.abxRisk) state.abxRisk[k].sort((a, b) => (b.prob || 0) - (a.prob || 0));
}

function openAbxModal(i) {
  const a = state.abById[i]; if (!a) return;
  const resist = state.abxResist[a.name] || [], risk = state.abxRisk[a.name] || [];
  const chips = [`<span class="abx-chip">AWaRe: <b style="color:${tierC(a.tier)}">${a.tier_name}</b></span>`];
  if (a.complexity != null) chips.push(`<span class="abx-chip">Complexity: <b>${a.complexity > 0 ? '+' : ''}${a.complexity.toFixed(2)}</b></span>`);
  if (a.auc != null) chips.push(`<span class="abx-chip">Predictability: <b>${a.auc.toFixed(2)} AUC</b></span>`);
  chips.push(`<span class="abx-chip">Resisted by <b>${resist.length}</b> species</span>`);
  const rlist = resist.slice(0, 10).map((r, k) =>
    `<div class="abx-li"><span class="rk">${k + 1}</span><span class="nm">${italic(r.p)}</span><span class="val res">${r.propR != null ? (r.propR * 100).toFixed(0) + '%' : '—'}</span></div>`).join('')
    || `<div class="abx-empty">No species are network&#8209;resistant to this drug yet.</div>`;
  const klist = risk.slice(0, 10).map((r, k) =>
    `<div class="abx-li"><span class="rk">${k + 1}</span><span class="nm">${italic(r.p)}</span><span class="val hot">${(r.prob * 100).toFixed(0)}%</span></div>`).join('')
    || `<div class="abx-empty">No at&#8209;risk species flagged for this drug.</div>`;
  $('#abxBody').innerHTML =
    `<div class="abx-head"><div class="abx-title">${a.name}</div></div>
     <div class="abx-chips">${chips.join('')}</div>
     <div class="abx-cols">
       <div class="abx-col"><h4>🔴 Most&#8209;resistant species</h4>${rlist}</div>
       <div class="abx-col"><h4>⚡ Most at risk next</h4>${klist}</div>
     </div>
     <p class="abx-note">"Resistant" = revealed comparative advantage ≥ 1 in the co&#8209;resistance network; the percentage is
     the latest observed share of resistant isolates. "At risk" = calibrated probability of acquiring resistance next.</p>`;
  $('#abxModal').hidden = false;
}

function renderEmptyCard() {
  const picks = ['Klebsiella pneumoniae', 'Escherichia coli', 'Acinetobacter baumannii',
    'Pseudomonas aeruginosa', 'Staphylococcus aureus', 'Neisseria gonorrhoeae'];
  const qp = picks.filter(p => state.detail[p]).map(p =>
    `<span class="qp" data-p="${p}">${italic(p)}</span>`).join('');
  $('#pathcard').className = 'pathcard empty';
  $('#pathcard').innerHTML = `<p class="hint">Start typing, or try a quick pick:</p><div class="quickpicks">${qp}</div>`;
}

// ============================================================ network
let VP; // viewport <g> for zoom/pan
const NODES = {}; // idx -> {g, core, halo, a, cx, cy}
function buildNetwork() {
  const svg = $('#net');
  svg.innerHTML = '';
  const W = 1000, H = 700, pad = 40;
  const innet = state.abx.map((a, i) => ({ a, i })).filter(o => o.a.x != null);
  const maxDeg = Math.max(...innet.map(o => o.a.degree || 0)) || 1;
  const X = x => pad + x * (W - 2 * pad);
  const Y = y => pad + (1 - y) * (H - 2 * pad);
  const R = d => 6 + Math.sqrt((d || 0) / maxDeg) * 12;

  // defs: soft glow
  const defs = el('defs');
  defs.innerHTML =
    `<filter id="glow" x="-80%" y="-80%" width="260%" height="260%">
       <feGaussianBlur stdDeviation="6" result="b"/><feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge>
     </filter>`;
  svg.appendChild(defs);

  VP = el('g', { id: 'viewport' });
  svg.appendChild(VP);

  // edges
  const eg = el('g');
  state.net.edges.forEach(([s, t, w]) => {
    const A = state.abById[s], B = state.abById[t];
    if (!A || A.x == null || !B || B.x == null) return;
    eg.appendChild(el('line', {
      class: 'edge', x1: X(A.x), y1: Y(A.y), x2: X(B.x), y2: Y(B.y),
      'stroke-width': 0.6 + w * 3.4, 'stroke-opacity': 0.18 + w * 0.5
    }));
  });
  VP.appendChild(eg);

  // nodes
  const ng = el('g');
  const built = [];
  innet.forEach(({ a, i }) => {
    const cx = X(a.x), cy = Y(a.y), r = R(a.degree);
    const halo = el('circle', { class: 'halo', cx, cy, r: r + 4, fill: '#ff3b2e', opacity: 0, filter: 'url(#glow)' });
    const core = el('circle', { class: 'core', cx, cy, r, fill: tierC(a.tier), stroke: '#0a1018', 'stroke-width': 1.4 });
    const label = el('text', { x: cx, y: cy - r - 5, 'text-anchor': 'middle', text: a.name });
    const g = el('g', { class: 'node', 'data-i': i }, [halo, core, label]);
    g.addEventListener('mouseenter', ev => showTip(i, ev));
    g.addEventListener('mousemove', ev => moveTip(ev));
    g.addEventListener('mouseleave', hideTip);
    g.addEventListener('click', ev => { ev.stopPropagation(); openAbxModal(i); });
    ng.appendChild(g);
    NODES[i] = { g, core, halo, a, cx, cy, r, label };
    built.push({ i, a, cx, cy, r, label });
  });
  // greedy label thinning: label the most-connected drugs first, drop any label whose
  // box would overlap one already kept. Dropped labels reappear on hover / when lit.
  const placed = [];
  built.sort((p, q) => (q.a.degree || 0) - (p.a.degree || 0));
  built.forEach(b => {
    const w = b.a.name.length * 5.9 + 6, h = 15;
    const box = { x: b.cx - w / 2, y: b.cy - b.r - 18, w, h };
    const hit = placed.some(o => box.x < o.x + o.w && box.x + box.w > o.x && box.y < o.y + o.h && box.y + box.h > o.y);
    if (hit) b.label.classList.add('lbl-off');
    else placed.push(box);
  });
  VP.appendChild(ng);

  initZoomPan(svg);
}

// spotlight a node in the network (triggered from the prediction list)
function focusNode(name) {
  clearFocus();
  const i = state.byName[name]; if (i == null) return;
  const N = NODES[i]; if (!N) return;
  N.g.classList.add('focus');
  N.g.parentNode.appendChild(N.g); // bring to front
  const ring = el('circle', { class: 'focusring', cx: N.cx, cy: N.cy, r: N.r + 7, fill: 'none', stroke: '#fff', 'stroke-width': 2.5 });
  ring.style.transformOrigin = `${N.cx}px ${N.cy}px`;
  VP.appendChild(ring);
  state._focusRing = ring; state._focusNode = N;
}
function clearFocus() {
  if (state._focusRing) { state._focusRing.remove(); state._focusRing = null; }
  if (state._focusNode) { state._focusNode.g.classList.remove('focus'); state._focusNode = null; }
}

let resetView = () => {};
function initZoomPan(svg) {
  let k = 1, tx = 0, ty = 0, dragging = false, sx, sy;
  const apply = () => VP.setAttribute('transform', `translate(${tx},${ty}) scale(${k})`);
  resetView = () => { k = 1; tx = 0; ty = 0; apply(); };
  svg.addEventListener('wheel', e => {
    e.preventDefault();
    const rect = svg.getBoundingClientRect();
    const mx = (e.clientX - rect.left) / rect.width * 1000, my = (e.clientY - rect.top) / rect.height * 700;
    // never zoom out past the initial default (k = 1); allow zooming in up to 5x
    const f = e.deltaY < 0 ? 1.12 : 1 / 1.12, nk = Math.max(1, Math.min(5, k * f));
    tx = mx - (mx - tx) * (nk / k); ty = my - (my - ty) * (nk / k); k = nk; apply();
  }, { passive: false });
  svg.addEventListener('mousedown', e => {
    if (e.target.closest('.node')) return; dragging = true; svg.classList.add('grabbing');
    sx = e.clientX; sy = e.clientY;
  });
  window.addEventListener('mousemove', e => {
    if (!dragging) return; const rect = svg.getBoundingClientRect();
    tx += (e.clientX - sx) / rect.width * 1000; ty += (e.clientY - sy) / rect.height * 700;
    sx = e.clientX; sy = e.clientY; apply();
  });
  window.addEventListener('mouseup', () => { dragging = false; svg.classList.remove('grabbing'); });
}

// ---------- tooltip ----------
function showTip(i, ev) {
  const a = state.abById[i], tip = $('#tooltip');
  let rows = `<div class="tt-row"><span>AWaRe tier</span><b style="color:${tierC(a.tier)}">${a.tier_name}</b></div>`;
  if (a.complexity != null) rows += `<div class="tt-row"><span>Complexity</span><b>${a.complexity > 0 ? '+' : ''}${a.complexity.toFixed(2)}</b></div>`;
  if (a.auc != null) rows += `<div class="tt-row"><span>Predictability (AUC)</span><b>${a.auc.toFixed(2)}</b></div>`;
  if (state.sel) {
    const d = state.detail[state.sel], ts = d.timeline[a.name];
    if (ts && ts.length) {
      const last = ts[ts.length - 1];
      const lit = d.lit[i] != null;
      rows += `<div class="tt-row"><span>${italic(state.sel)}</span><b class="${lit ? 'tt-hot' : ''}">${(last[1] * 100).toFixed(0)}% resistant (${last[0]})</b></div>`;
      if (lit) rows += `<div class="tt-row"><span>Network-resistant since</span><b class="tt-hot">${d.lit[i]}</b></div>`;
    } else {
      rows += `<div class="tt-row"><span>${italic(state.sel)}</span><b style="color:#6d8299">not tested</b></div>`;
    }
    highlightHeatRow(a.name);
  }
  tip.innerHTML = `<div class="tt-name">${a.name}</div>${rows}`;
  tip.hidden = false; moveTip(ev);
}
function moveTip(ev) {
  const stage = $('.net-stage').getBoundingClientRect(), tip = $('#tooltip');
  const node = ev.currentTarget.querySelector('.core').getBoundingClientRect();
  tip.style.left = (node.left + node.width / 2 - stage.left) + 'px';
  tip.style.top = (node.top - stage.top) + 'px';
}
function hideTip() { $('#tooltip').hidden = true; clearHeatHl(); }

// ============================================================ search
function initSearch() {
  const inp = $('#search'), sug = $('#suggest');
  let hi = -1, list = [];
  const render = q => {
    q = q.trim().toLowerCase();
    list = state.paths
      .filter(p => p.name.toLowerCase().includes(q))
      .sort((a, b) => (b.n_records - a.n_records));
    const items = list.map((p, i) =>
      `<div class="item${i === hi ? ' hi' : ''}" data-p="${p.name}">
         <span class="nm">${italic(p.name)}</span>
         ${p.esk ? `<span class="tag esk">${p.esk}</span>` : ''}
       </div>`).join('');
    const count = `<div class="sug-count">${list.length} of ${state.paths.length} species</div>`;
    sug.innerHTML = list.length ? count + items : `<div class="item"><span class="nm" style="color:#6d8299">No match</span></div>`;
    sug.classList.add('on');
  };
  inp.addEventListener('input', () => { hi = -1; render(inp.value); });
  inp.addEventListener('focus', () => render(inp.value));
  inp.addEventListener('keydown', e => {
    if (e.key === 'ArrowDown') { hi = Math.min(hi + 1, list.length - 1); render(inp.value); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { hi = Math.max(hi - 1, 0); render(inp.value); e.preventDefault(); }
    else if (e.key === 'Enter' && list[hi < 0 ? 0 : hi]) { pick(list[hi < 0 ? 0 : hi].name); }
    else if (e.key === 'Escape') { sug.classList.remove('on'); }
  });
  sug.addEventListener('mousedown', e => { const it = e.target.closest('.item'); if (it && it.dataset.p) pick(it.dataset.p); });
  document.addEventListener('click', e => { if (!e.target.closest('.search-wrap')) sug.classList.remove('on'); });
  function pick(name) { inp.value = ''; sug.classList.remove('on'); inp.blur(); selectPathogen(name); }
}

// ============================================================ select pathogen
function selectPathogen(name) {
  if (!state.detail[name]) return;
  state.sel = name;
  clearFocus();
  try { history.replaceState(null, '', `?p=${encodeURIComponent(name)}`); } catch (e) {}
  const p = state.paths.find(x => x.name === name), d = state.detail[name];
  renderCard(p, d);
  renderForecast(p, d);
  renderHeat(name, d);
  // network -> lighting mode
  $('#net-title').innerHTML = `How <i>${name}</i> conquers the network`;
  $('#net-sub').textContent = 'Nodes ignite the year this pathogen becomes network-resistant to that drug. Press play to watch it spread.';
  $('#timeplay').hidden = false;
  $('#spreadbadge').hidden = false;
  $('#deselect').hidden = false;
  const firstYears = Object.values(d.lit);
  const start = firstYears.length ? Math.min(...firstYears) : state.meta.year_min;
  d._start = start;
  setYear(state.meta.year_max);
  $('#yearslider').value = state.meta.years.indexOf(state.meta.year_max);
}

// return to the initial AWaRe-coloured network
function deselectPathogen() {
  stopPlay();
  state.sel = null;
  try { history.replaceState(null, '', location.pathname); } catch (e) {}
  hideTip();
  // hide the pathogen-only UI
  $('#timeplay').hidden = true;
  $('#spreadbadge').hidden = true;
  $('#deselect').hidden = true;
  closeAllModals();
  renderEmptyCard();
  $('#forecast').innerHTML = '';
  // restore network chrome + node visuals + view
  $('#net-title').innerHTML = 'Antibiotic co&#8209;resistance network';
  $('#net-sub').textContent = 'Each node is an antibiotic; links join drugs that tend to fail together. Colour = WHO AWaRe tier.';
  clearFocus();
  resetNetworkVisual();
  resetView();
}

// clear per-pathogen lighting; back to tier colours, no glow, no dimming
function resetNetworkVisual() {
  for (const i in NODES) {
    const N = NODES[i];
    N.g.classList.remove('lit', 'dim', 'sel');
    N.halo.setAttribute('opacity', 0);
    N.halo.setAttribute('r', (N.r + 4).toFixed(1));
    N.core.setAttribute('fill', tierC(N.a.tier));
  }
}

function renderCard(p, d) {
  const badges = [];
  if (p.esk) badges.push(`<span class="badge esk">${p.esk} priority</span>`);
  const priLab = { 3: 'WHO Critical', 2: 'WHO High', 1: 'WHO Medium' }[p.priority];
  if (priLab) badges.push(`<span class="badge pri${p.priority}">${priLab}</span>`);
  if (!badges.length) badges.push(`<span class="badge neutral">Environmental / non-priority</span>`);
  const fitPct = p.fitness_pct != null ? Math.round(p.fitness_pct) : null;
  const nLit = Object.keys(d.lit).length;
  $('#pathcard').className = 'pathcard';
  $('#pathcard').innerHTML = `
    <div class="pc-name">${p.name}</div>
    <div class="pc-badges">${badges.join('')}</div>
    <div class="pc-stats">
      <div class="pcs"><div class="v">${nLit}</div><div class="k">drugs resisted (network)</div></div>
      <div class="pcs"><div class="v">${p.n_abx}</div><div class="k">drugs ever tested</div></div>
    </div>
    ${fitPct != null ? `<div class="fitbar-wrap">
      <div class="fitbar-lab"><span>Resistance fitness</span><span><b style="color:#eaf1f8">${fitPct}ᵗʰ</b> percentile</span></div>
      <div class="fitbar"><div style="width:${fitPct}%"></div></div>
    </div>` : ''}
    ${Object.keys(d.timeline).length ? `<button class="timelinebtn" id="openHeat">📈 View resistance timeline</button>` : ''}`;
}

function renderForecast(p, d) {
  const box = $('#forecast');
  const fc = (d.forecast || []).filter(f => f.prob != null);
  if (!fc.length) {
    box.innerHTML = `<h3>⚡ Next-move forecast</h3><p class="fc-empty">No susceptible drugs with enough network signal to forecast for this pathogen.</p>`;
    return;
  }
  const maxP = Math.max(...fc.map(f => f.prob));
  const rows = fc.slice(0, 6).map((f, i) => {
    const ai = state.byName[f.abx], a = ai != null ? state.abById[ai] : null;
    const tier = a ? a.tier_name : '';
    return `<div class="fc-row" data-abx="${f.abx}">
      <div class="fc-rank">${i + 1}</div>
      <div class="fc-main"><div class="fc-name">${f.abx}</div><div class="fc-tier">${tier}${f.cur_propR != null ? ` · now ${(f.cur_propR * 100).toFixed(0)}% resistant` : ''}</div></div>
      <div class="fc-meter"><div class="fc-bar"><div style="width:${Math.max(6, f.prob / maxP * 100)}%"></div></div>
        <div class="fc-pct">${(f.prob * 100).toFixed(0)}%</div></div>
    </div>`;
  }).join('');
  box.innerHTML = `<h3>⚡ Next-move forecast</h3>
    <p class="sub">Estimated chance <i>${p.name}</i> acquires resistance to each drug it can still be treated with, from its position in the network. <b style="color:#9fb2c6">Click a drug to locate it on the network.</b></p>${rows}`;
  box.querySelectorAll('.fc-row').forEach(r => {
    r.addEventListener('click', () => focusNode(r.dataset.abx));
  });
}

// ============================================================ time / animation
function initTime() {
  const sl = $('#yearslider');
  sl.max = state.meta.years.length - 1;
  sl.addEventListener('input', () => { stopPlay(); setYear(state.meta.years[+sl.value]); });
  $('#playbtn').addEventListener('click', togglePlay);
}
function setYear(y) {
  state.year = y;
  $('#yearval').textContent = y;
  const idx = state.meta.years.indexOf(y);
  if (+$('#yearslider').value !== idx) $('#yearslider').value = idx;
  if (!state.sel) return;
  const d = state.detail[state.sel];
  let count = 0;
  for (const i in NODES) {
    const litYear = d.lit[i];
    const lit = litYear != null && litYear <= y;
    const N = NODES[i];
    N.g.classList.toggle('lit', lit);
    N.g.classList.toggle('dim', !lit);
    // intensity by this pathogen's resistance fraction near this year
    let inten = 0.5;
    const ts = d.timeline[N.a.name];
    if (ts) { const rec = lastAtOrBefore(ts, y); if (rec) inten = rec[1]; }
    N.halo.setAttribute('opacity', lit ? (0.35 + inten * 0.5).toFixed(2) : 0);
    N.halo.setAttribute('r', (N.r + 4 + (lit ? inten * 7 : 0)).toFixed(1));
    N.core.setAttribute('fill', lit ? tierC(N.a.tier) : '#33465c');
    if (lit) count++;
  }
  const badge = $('#spreadbadge');
  badge.innerHTML = `🔴 ${count} / ${Object.keys(NODES).length} drugs defeated by ${state.year}`;
}
function lastAtOrBefore(ts, y) { let out = null; for (const r of ts) { if (r[0] <= y) out = r; else break; } return out; }

function togglePlay() { state.playing ? stopPlay() : startPlay(); }
function startPlay() {
  if (!state.sel) return;
  const d = state.detail[state.sel], years = state.meta.years;
  let idx = years.indexOf(state.year);
  if (idx >= years.length - 1) { idx = Math.max(0, years.indexOf(d._start) - 1); } // restart from just before first ignition
  $('#playbtn').textContent = '❚❚';
  state.playing = setInterval(() => {
    idx++;
    if (idx >= years.length) { stopPlay(); return; }
    setYear(years[idx]);
  }, 620);
}
function stopPlay() { if (state.playing) { clearInterval(state.playing); state.playing = null; } $('#playbtn').textContent = '▶'; }

// ============================================================ heatmap
let HEAT_ROWS = {};
function renderHeat(name, d) {
  const svg = $('#heat'); svg.innerHTML = ''; HEAT_ROWS = {};
  $('#heat-name').innerHTML = italic(name);
  const years = state.meta.years;
  const abxNames = Object.keys(d.timeline);
  if (!abxNames.length) return;
  // order: network-resistant first (by ignition year), then the rest alphabetically
  const litYear = a => { const i = state.byName[a]; return i != null && d.lit[i] != null ? d.lit[i] : Infinity; };
  abxNames.sort((a, b) => (litYear(a) - litYear(b)) || a.localeCompare(b));

  const padL = 168, padT = 8, padB = 26, cw = 30, ch = 19;
  const W = padL + years.length * cw + 8, Hh = padT + abxNames.length * ch + padB;
  svg.setAttribute('viewBox', `0 0 ${W} ${Hh}`);
  svg.setAttribute('width', W); svg.setAttribute('height', Hh);
  const g = el('g');

  // year labels
  years.forEach((y, j) => {
    if (y % 3 === 1 || j === years.length - 1) {
      g.appendChild(el('text', { class: 'hm-col-lab', x: padL + j * cw + cw / 2, y: Hh - 10, 'text-anchor': 'middle', text: y }));
    }
  });

  abxNames.forEach((a, r) => {
    const yPix = padT + r * ch;
    const ai = state.byName[a];
    const lit = ai != null && d.lit[ai] != null;
    g.appendChild(el('text', {
      class: 'hm-row-lab', x: padL - 8, y: yPix + ch / 2 + 3, 'text-anchor': 'end',
      text: a.length > 24 ? a.slice(0, 23) + '…' : a, fill: lit ? '#ffb9a3' : '#8ea3ba'
    }));
    const rowCells = [];
    years.forEach((y, j) => {
      const rec = d.timeline[a].find(t => t[0] === y);
      const fill = rec ? magma(rec[1]) : '#152234';
      const cell = el('rect', { class: 'hm-cell', x: padL + j * cw, y: yPix, width: cw - 1.5, height: ch - 1.5, fill, rx: 2 });
      if (rec) { cell.setAttribute('data-tip', `${(rec[1] * 100).toFixed(0)}% of ${rec[2]} tested (${y})`); }
      rowCells.push(cell);
      g.appendChild(cell);
    });
    HEAT_ROWS[a] = rowCells;
  });
  svg.appendChild(g);
  $('#heat-legend').innerHTML = `<span>0%</span><span class="grad"></span><span>100% resistant</span><span style="margin-left:14px"><span style="display:inline-block;width:11px;height:11px;background:#152234;border-radius:2px;vertical-align:-1px"></span> no data</span>`;
}
function highlightHeatRow(a) {
  clearHeatHl();
  (HEAT_ROWS[a] || []).forEach(c => c.classList.add('hl'));
}
function clearHeatHl() { document.querySelectorAll('.hm-cell.hl').forEach(c => c.classList.remove('hl')); }
