// js/components/reaction.js
// RDKit.js-based SMIRKS reaction rendering for EVODEX.2.
//
// Exports:
//   window.renderReactionInto(container, spec)
//
// Where `spec` can be:
//   - a string SMIRKS, or
//   - an object with `smirks` (or `reaction_smirks`) field.
//
// This renders the full reaction in a single RDKit-generated SVG.

(function (global) {
    "use strict";
  
    let RDKit = null;
    let rdkitReady = null;
  
    // ---------------------------------------------------------------------------
    // Lazy RDKit initialization
    // ---------------------------------------------------------------------------
  
    function ensureRDKit() {
      if (RDKit) return Promise.resolve(RDKit);
      if (rdkitReady) return rdkitReady;
  
      if (typeof initRDKitModule !== "function") {
        console.warn("RDKit_minimal.js not yet loaded: initRDKitModule is undefined.");
        return Promise.resolve(null);
      }
  
      rdkitReady = initRDKitModule()
        .then((instance) => {
          RDKit = instance;
          return instance;
        })
        .catch((err) => {
          console.error("Failed to initialize RDKit:", err);
          RDKit = null;
          return null;
        });
  
      return rdkitReady;
    }
  
    // ---------------------------------------------------------------------------
    // Helper: extract a SMIRKS string from whatever the caller passed
    // ---------------------------------------------------------------------------
  
    function extractSmirks(spec) {
      if (!spec) return null;
  
      if (typeof spec === "string") {
        const s = spec.trim();
        return s.length ? s : null;
      }
  
      if (typeof spec === "object") {
        if (spec.smirks && String(spec.smirks).trim()) {
          return String(spec.smirks).trim();
        }
        if (spec.reaction_smirks && String(spec.reaction_smirks).trim()) {
          return String(spec.reaction_smirks).trim();
        }
      }
  
      return null;
    }
  
    // ---------------------------------------------------------------------------
    // Main entry: render a reaction SMIRKS into `container`
    // ---------------------------------------------------------------------------
  
    function renderReactionInto(container, spec) {
      if (!container) return;
  
      container.innerHTML = "";
  
      const smirks = extractSmirks(spec);
      if (!smirks) {
        container.textContent = "Reaction SMIRKS not available.";
        return;
      }
  
      // Placeholder while RDKit / WASM gets ready
      const placeholder = document.createElement("div");
      placeholder.textContent = "Rendering reaction…";
      container.appendChild(placeholder);
  
      ensureRDKit().then((instance) => {
        if (!instance) {
          container.innerHTML = "";
          const pre = document.createElement("pre");
          pre.textContent = smirks;
          container.appendChild(
            document.createTextNode("RDKit is not available; showing SMIRKS:")
          );
          container.appendChild(document.createElement("br"));
          container.appendChild(pre);
          return;
        }
  
        container.innerHTML = "";
  
        try {
          const rxn = instance.get_rxn(smirks);
          if (!rxn) {
            container.textContent = "Failed to create RDKit reaction from SMIRKS.";
            return;
          }
  
          // Some builds expose is_valid(); if present, use it
          if (typeof rxn.is_valid === "function" && !rxn.is_valid()) {
            rxn.delete();
            container.textContent = "Invalid SMIRKS.";
            return;
          }
  
          // Choose width based on container; fallback if necessary
          const padding = 40;
          const cw = container.clientWidth || 700;
          const width = Math.max(300, cw - padding);
          const height = 220;
  
          const svgText = rxn.get_svg(width, height);
          rxn.delete();
  
          // Insert SVG
          container.innerHTML = svgText;
  
          // Make the SVG responsive (max-width: 100%)
          const svg = container.querySelector("svg");
          if (svg) {
            svg.removeAttribute("width");
            svg.removeAttribute("height");
            svg.setAttribute("width", "100%");
            svg.setAttribute("height", "auto");
            svg.style.display = "block";
          }
        } catch (err) {
          console.error("Error rendering RDKit reaction from SMIRKS:", err);
          container.textContent = "Error rendering reaction.";
        }
      });
    }
  
    // Expose globally for operator_detail.js (and others)
    global.renderReactionInto = renderReactionInto;
  })(window);