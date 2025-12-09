// js/router.js
// Simple hash router for the EVODEX.2 single page app.
//
// Routes:
//
//   #/                      -> home
//   #/type/E                -> type index for EVODEX-E
//   #/type/P                -> type index for EVODEX-P
//   #/E/EVODEX.2-E1         -> operator detail page
//   #/P/EVODEX.2-P1240      -> operator detail page
//
// Unknown routes fall back to the home view.

/**
 * Parse the current location.hash into a route description.
 *
 * Returns an object:
 *   { view: "home" }
 *   { view: "typeIndex", type: "E" }
 *   { view: "operatorDetail", type: "E", id: "EVODEX.2-E1" }
 */
function parseHash() {
    let hash = window.location.hash || "";
    hash = hash.replace(/^#/, "").trim();  // remove leading '#'
  
    if (!hash || hash === "/") {
      return { view: "home" };
    }
  
    // Remove leading slash if present
    if (hash.startsWith("/")) {
      hash = hash.slice(1);
    }
  
    const parts = hash.split("/");
  
    // #/type/E
    if (parts[0] === "type" && parts[1]) {
      return { view: "typeIndex", type: parts[1] };
    }
  
    // #/E/EVODEX.2-E1
    if (parts.length === 2 && parts[0] && parts[1]) {
      return { view: "operatorDetail", type: parts[0], id: parts[1] };
    }
  
    // Fallback
    return { view: "home" };
  }
  
  /**
   * Render the current route into the #app element.
   */
  function renderRoute() {
    const app = document.getElementById("app");
    if (!app) return;
  
    const route = parseHash();
  
    if (route.view === "home") {
      if (typeof renderHome === "function") {
        renderHome(app);
      } else {
        app.innerHTML = `
          <section>
            <h2>EVODEX.2</h2>
            <p>Home view is not implemented yet.</p>
          </section>
        `;
      }
      return;
    }
  
    if (route.view === "typeIndex") {
      const type = route.type;
      if (typeof renderTypeIndex === "function") {
        renderTypeIndex(app, type);
      } else {
        app.innerHTML = `
          <section>
            <h2>EVODEX-${type}</h2>
            <p>Type index view is not implemented yet.</p>
          </section>
        `;
      }
      return;
    }
  
    if (route.view === "operatorDetail") {
      const type = route.type;
      const id = route.id;
      if (typeof renderOperatorDetail === "function") {
        renderOperatorDetail(app, type, id);
      } else {
        app.innerHTML = `
          <section>
            <h2>${id}</h2>
            <p>Operator detail view is not implemented yet.</p>
          </section>
        `;
      }
      return;
    }
  
    // Final fallback
    app.innerHTML = `
      <section>
        <h2>Not found</h2>
        <p>The requested view could not be resolved.</p>
      </section>
    `;
  }
  
  /**
   * Initialize the app once the DOM is ready.
   * Renders the header and the initial route, then attaches hashchange listener.
   */
  function initEVODEXApp() {
    const headerEl = document.getElementById("site-header");
    if (headerEl && typeof renderHeader === "function") {
      renderHeader(headerEl);
    }
  
    renderRoute();
    window.addEventListener("hashchange", renderRoute);
  }
  
  window.addEventListener("DOMContentLoaded", initEVODEXApp);