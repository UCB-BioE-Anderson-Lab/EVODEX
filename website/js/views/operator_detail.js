// js/views/operator_detail.js
// Operator detail view, e.g. #/E/EVODEX.2-E1

function renderOperatorDetail(app, rawType, id) {
  const type = normalizeTypeForDetail(rawType);

  if (!type) {
    app.innerHTML = card("Operator detail", `<p>Unknown EVODEX type: ${rawType}</p>`);
    return;
  }

  if (!id) {
    app.innerHTML = card("Operator detail", `<p>No operator ID specified.</p>`);
    return;
  }

  if (typeof EVODEXData === "undefined" || !EVODEXData.TYPE_CONFIG[type]) {
    app.innerHTML = card(
      "Operator detail",
      `<p>No data configuration found for type: ${escapeHtml(type)}.</p>`
    );
    return;
  }

  app.innerHTML = card(id, `<p>Loading EVODEX-${escapeHtml(type)} entry…</p>`);

  EVODEXData.getEntry(type, id)
    .then(function (entry) {
      if (!entry) {
        app.innerHTML = card(id, `<p>No data found for <code>${escapeHtml(id)}</code>.</p>`);
        return;
      }

      // Build the page using DOM nodes so we can append multiple cards.
      app.innerHTML = "";
      var frag = document.createDocumentFragment();

      // Header card: ID + type
      frag.appendChild(
        cardNode(id, `
          <p><strong>Type:</strong> EVODEX-${escapeHtml(type)}</p>
        `)
      );

      // Reaction / SMIRKS card (SMIRKS-centric)
      if (hasSmirks(entry)) {
        var reactionSection = buildReactionSection(entry);
        frag.appendChild(reactionSection);
      }

      // F-type: formula operators
      if (type === "F") {
        var fSection = buildFormulaSection(entry);
        if (fSection) frag.appendChild(fSection);
      }

      // M-type: mass operators
      if (type === "M" || type === "M_SUBSET") {
        var mSection = buildMassSection(entry);
        if (mSection) frag.appendChild(mSection);
      }

      // Sources card (for any type that has a sources field)
      var sourcesSection = buildSourcesSection(entry);
      if (sourcesSection) frag.appendChild(sourcesSection);

      app.appendChild(frag);
    })
    .catch(function (err) {
      console.error(err);
      app.innerHTML = card(
        id,
        `<p>Failed to load operator data for <code>${escapeHtml(id)}</code>.</p>`
      );
    });
}

/* ---------- Helpers ---------- */

// Roughly mirror normalizeTypeForView so we can accept "e", "E", "am", "Am", etc.
function normalizeTypeForDetail(rawType) {
  if (!rawType) return null;
  var t = String(rawType).trim();

  if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[t]) return t;

  var upper = t.toUpperCase();
  if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[upper]) return upper;

  var cap = t.charAt(0).toUpperCase() + t.slice(1);
  if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[cap]) return cap;

  return null;
}

function hasSmirks(entry) {
  return entry && typeof entry.smirks === "string" && entry.smirks.trim().length > 0;
}

// Kept for now (unused with SMIRKS-based rendering, but harmless)
function hasReactionSides(entry) {
  return entry &&
    (Array.isArray(entry.substrates) || Array.isArray(entry.products));
}

// Build the SMIRKS + reaction card
function buildReactionSection(entry) {
  var container = cardNode("Reaction", "");
  var body = container.querySelector(".evodex-card-body") || container;

  // Reaction SVG via RDKit (from full SMIRKS)
  if (hasSmirks(entry) && typeof renderReactionInto === "function") {
    var reactionDiv = document.createElement("div");
    reactionDiv.className = "evodex-reaction-container";
    body.appendChild(reactionDiv);

    // New semantics: pass full SMIRKS to reaction.js
    renderReactionInto(reactionDiv, entry.smirks);
  }

  // SMIRKS text
  if (hasSmirks(entry)) {
    var smirksBlock = document.createElement("pre");
    smirksBlock.className = "smirks-text";
    smirksBlock.textContent = entry.smirks;
    body.appendChild(smirksBlock);
  }

  return container;
}

// Build F-type (formula) card
function buildFormulaSection(entry) {
  if (!entry) return null;

  var formulaDiff = entry.formula_diff || entry.formulaDiff || entry.formula;
  if (!formulaDiff) return null;

  var html = `
    <p><strong>Formula difference</strong></p>
    <pre class="smirks-text">${escapeHtml(String(formulaDiff))}</pre>
  `;
  return cardNode("Formula operator (F)", html);
}

// Build M-type (mass) card
function buildMassSection(entry) {
  if (!entry) return null;

  var mass = entry.mass;
  if (mass == null) return null;

  var html = `
    <p><strong>Exact mass difference (Da)</strong></p>
    <p>${escapeHtml(String(mass))}</p>
  `;

  if (entry.formula) {
    html += `
      <p><strong>Associated formula</strong></p>
      <pre class="smirks-text">${escapeHtml(String(entry.formula))}</pre>
    `;
  }

  return cardNode("Mass operator (M)", html);
}

// Build Sources card (if any)
function buildSourcesSection(entry) {
  if (!entry || entry.sources == null) return null;

  var sourceIds = parseSources(entry.sources);
  if (!sourceIds.length) return null;

  var itemsHtml = sourceIds.map(function (src) {
    var srcId = String(src).trim();
    if (!srcId) return "";

    var srcType = inferTypeFromId(srcId) || "";
    var typeSegment = srcType || "";

    if (typeSegment) {
      var href = "#/" + typeSegment + "/" + srcId;
      return (
        '<li><a href="' +
        href +
        '">' +
        escapeHtml(srcId) +
        "</a></li>"
      );
    } else {
      return "<li>" + escapeHtml(srcId) + "</li>";
    }
  }).join("");

  var html = "<ul class=\"evodex-id-list\">" + itemsHtml + "</ul>";
  return cardNode("Sources", html);
}

// sources may be array or comma-separated string
function parseSources(raw) {
  if (!raw) return [];
  if (Array.isArray(raw)) return raw.slice();
  if (typeof raw === "string") {
    return raw.split(",").map(function (s) { return s.trim(); }).filter(Boolean);
  }
  try {
    var parsed = JSON.parse(raw);
    if (Array.isArray(parsed)) return parsed;
  } catch (e) {
    // fall through
  }
  return [];
}

// Try to infer the operator type (E, P, R, F, M, ...) from an ID like "EVODEX.2-E1"
function inferTypeFromId(id) {
  var m = /EVODEX\.2-([A-Za-z]+)/.exec(id);
  if (!m) return null;
  var t = m[1];
  if (EVODEXData && EVODEXData.TYPE_CONFIG) {
    if (EVODEXData.TYPE_CONFIG[t]) return t;
    var upper = t.toUpperCase();
    if (EVODEXData.TYPE_CONFIG[upper]) return upper;
    var cap = t.charAt(0).toUpperCase() + t.slice(1);
    if (EVODEXData.TYPE_CONFIG[cap]) return cap;
  }
  return t;
}

// Minimal HTML escaper for safety
function escapeHtml(str) {
  return String(str)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");
}