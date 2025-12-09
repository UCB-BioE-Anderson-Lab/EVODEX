// js/components/header.js
// Renders the EVODEX header bar, matching the EVODEX.1 style,
// but cleaned slightly for responsiveness and SPA behavior.
//

function renderHeader(el) {
    if (!el) return;
  
    el.innerHTML = `
      <div class="evodex-header-inner">
        <a href="#/" class="evodex-logo-link">
          <img src="website/images/evodex_logo.png" class="evodex-logo" alt="EVODEX logo">
        </a>
        <div class="evodex-title-block">
          <h1 class="evodex-title">EVODEX.2</h1>
          <div class="evodex-subtitle">Evolutionary Decomposition of Chemical Transformations</div>
        </div>
      </div>
    `;
  }