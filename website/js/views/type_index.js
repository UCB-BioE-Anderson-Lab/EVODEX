// js/views/type_index.js
// Type index view (e.g. #/type/E)

function renderTypeIndex(app, rawType) {
  const type = normalizeTypeForView(rawType);

  if (!type) {
    app.innerHTML = card("Type Index", `<p>Unknown EVODEX type: ${rawType}</p>`);
    return;
  }

  if (typeof EVODEXData === "undefined" || !EVODEXData.TYPE_CONFIG[type]) {
    app.innerHTML = card("Type Index", `<p>No data configuration found for type: ${type}</p>`);
    return;
  }

  app.innerHTML = card(`EVODEX-${type} index`, `<p>Loading EVODEX-${type} entries…</p>`);

  EVODEXData.listIds(type)
    .then((ids) => {
      if (!ids || ids.length === 0) {
        app.innerHTML = card(
          `EVODEX-${type} index`,
          `<p>No entries found for EVODEX-${type}.</p>`
        );
        return;
      }

      // Sort by numeric suffix (last number in the ID), then by full string
      const sorted = ids.slice().sort(compareEvodexIds);

      const listHtml = `
        <ul class="evodex-id-list">
          ${sorted.map((id) => idListItem(type, id)).join("")}
        </ul>
      `;

      app.innerHTML = card(`EVODEX-${type} index`, listHtml);
    })
    .catch((err) => {
      console.error(err);
      app.innerHTML = card(
        `EVODEX-${type} index`,
        `<p>Failed to load EVODEX-${type} data.</p>`
      );
    });
}

// Normalize the type string for display / data lookup
function normalizeTypeForView(rawType) {
  if (!rawType) return null;
  const t = String(rawType).trim();

  if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[t]) return t;

  const upper = t.toUpperCase();
  if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[upper]) return upper;

  const cap = t.charAt(0).toUpperCase() + t.slice(1);
  if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[cap]) return cap;

  return null;
}

function idListItem(type, id) {
  const hashType = type; // already normalized
  return `
    <li>
      <a href="#/${hashType}/${id}">${id}</a>
    </li>
  `;
}

/**
 * Compare two EVODEX IDs in a natural "…E1, E2, … E10, E11…" order.
 * Uses the *last* numeric chunk in the string as the primary sort key.
 */
function compareEvodexIds(a, b) {
  const na = extractLastNumber(a);
  const nb = extractLastNumber(b);

  if (na != null && nb != null && na !== nb) {
    return na - nb;
  }

  // If no numbers or equal numerically, fall back to plain string compare
  return String(a).localeCompare(String(b));
}

/**
 * Extract the last integer found in a string, or null if none.
 * For "EVODEX.2-E123" this returns 123 (not 2).
 */
function extractLastNumber(id) {
  // (\d+)(?!.*\d) = "a run of digits not followed by any more digits later"
  const m = /(\d+)(?!.*\d)/.exec(String(id));
  if (!m) return null;
  const n = parseInt(m[1], 10);
  return Number.isNaN(n) ? null : n;
}