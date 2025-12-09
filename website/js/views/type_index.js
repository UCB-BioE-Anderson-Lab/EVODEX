// js/views/type_index.js
// Type index view (e.g. #/type/E)

function renderTypeIndex(app, rawType) {
    const type = normalizeTypeForView(rawType);
  
    if (!type) {
      app.innerHTML = card("Type Index", `<p>Unknown EVODEX type: ${rawType}</p>`);
      return;
    }
  
    // If EVODEXData or TYPE_CONFIG is missing, bail gracefully
    if (typeof EVODEXData === "undefined" || !EVODEXData.TYPE_CONFIG[type]) {
      app.innerHTML = card("Type Index", `<p>No data configuration found for type: ${type}</p>`);
      return;
    }
  
    // Initial loading state
    app.innerHTML = card(`EVODEX-${type} index`, `<p>Loading EVODEX-${type} entries…</p>`);
  
    EVODEXData.listIds(type)
      .then((ids) => {
        if (!ids || ids.length === 0) {
          app.innerHTML = card(`EVODEX-${type} index`, `<p>No entries found for EVODEX-${type}.</p>`);
          return;
        }
  
        const listHtml = `
          <ul class="evodex-id-list">
            ${ids.map((id) => idListItem(type, id)).join("")}
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
  
    // Try exact first (for Am, Cm, etc.)
    if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[t]) return t;
  
    // Then uppercase (E, P, C, etc.)
    const upper = t.toUpperCase();
    if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[upper]) return upper;
  
    // Capitalize (Am, Cm, etc.)
    const cap = t.charAt(0).toUpperCase() + t.slice(1);
    if (EVODEXData && EVODEXData.TYPE_CONFIG && EVODEXData.TYPE_CONFIG[cap]) return cap;
  
    return null;
  }
  
  function idListItem(type, id) {
    // Link format: #/E/EVODEX.2-E1, #/P/EVODEX.2-P123, etc.
    // The type in the hash is the "type" segment (E, P, A, Am, etc.).
    const hashType = type; // already normalized
    return `
      <li>
        <a href="#/${hashType}/${id}">${id}</a>
      </li>
    `;
  }