// js/views/home.js
// Home view for EVODEX.2

function renderHome(app) {
  // defensive: if card() isn't defined yet, fall back to very simple markup
  if (typeof card !== "function") {
    app.innerHTML = `
      <section>
        <h2>Home</h2>
        <p>EVODEX.2 single-page app is wired up, but the card component is missing.</p>
      </section>
    `;
    return;
  }

  const introHtml = `
    <p>
      Welcome to <strong>EVODEX.2</strong>, a dynamic viewer for EVODEX operators.
      Use the links below to explore operators by type.
    </p>
    <p>
      Click on a type to see its index, then click an entry ID to view the full
      operator page with dynamically rendered reactions.
    </p>
    <ul class="evodex-type-list">
      ${buildTypeListItems()}
    </ul>
  `;

  app.innerHTML = card("Home", introHtml);
}

// Build the list of EVODEX types to show on the home page.
function buildTypeListItems() {
  // Curated, ordered list of types we want to expose in the UI.
  // (We intentionally do NOT show helper types like E_OPS, SELECTED here.)
  const orderedTypes = [
    "R",    // full reactions (raw corpus)
    "P",    // partial reactions
    "E",    // E operators
    "A",
    "Am",
    "B",
    "Bm",
    "C",
    "Cm",
    "D",
    "Dm",
    "F",    // formula operators
    "M",    // mass operators
    "M_SUBSET"
  ];

  // TYPE_CONFIG comes from data_loader.js
  if (typeof EVODEXData === "undefined" || !EVODEXData.TYPE_CONFIG) {
    // fall back to static list if for some reason TYPE_CONFIG isn't ready
    return orderedTypes.map(t => typeListItem(t)).join("");
  }

  return orderedTypes
    .filter(t => EVODEXData.TYPE_CONFIG[t])   // only show types that really exist
    .map(t => typeListItem(t))
    .join("");
}

function typeListItem(type) {
  const labelMap = {
    R: "EVODEX-R (full reactions)",
    P: "EVODEX-P (partial reactions)",
    E: "EVODEX-E (operators)",
    A: "EVODEX-A",
    Am: "EVODEX-Am",
    B: "EVODEX-B",
    Bm: "EVODEX-Bm",
    C: "EVODEX-C",
    Cm: "EVODEX-Cm",
    D: "EVODEX-D",
    Dm: "EVODEX-Dm",
    F: "EVODEX-F (formula operators)",
    M: "EVODEX-M (mass operators)",
    M_SUBSET: "EVODEX-M (MS subset)"
  };
  const text = labelMap[type] || `EVODEX-${type}`;
  return `
    <li>
      <a href="#/type/${type}">${text}</a>
    </li>
  `;
}