// js/components/card.js
// A simple utility for generating EVODEX-style card sections.
//
// Usage:
//   app.innerHTML = card("Title", "<p>content</p>");
//   app.appendChild(cardNode("Title", "<p>content</p>"));
//
// This wrapper keeps markup consistent across all views.

function card(title, innerHTML) {
    const safeTitle = title ? `<h2>${title}</h2>` : "";
    return `
      <section class="evodex-card">
        ${safeTitle}
        <div class="evodex-card-body">
          ${innerHTML || ""}
        </div>
      </section>
    `;
  }
  
  // DOM-producing version — helpful when building multiple cards
  function cardNode(title, innerHTML) {
    const div = document.createElement("div");
    div.innerHTML = card(title, innerHTML).trim();
    return div.firstElementChild;
  }