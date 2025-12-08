// website/reaction_renderer.js
(function (global) {
    "use strict";
  
    if (typeof SmilesDrawer === "undefined") {
      console.error("SmilesDrawer must be loaded before reaction_renderer.js");
      return;
    }
  
    const drawerOptions = {
      width: 220,
      height: 140,
      padding: 15,
      explicitHydrogens: true,
      terminalCarbons: true,
      isometric: true
    };
  
    const drawer = new SmilesDrawer.Drawer(drawerOptions);

    function annotateIsotopeCarbons() {
      if (!drawer || !drawer.graph || !drawer.canvasWrapper || !drawer.canvasWrapper.canvas) {
        console.warn('annotateIsotopeCarbons: drawer/graph/canvasWrapper/canvas missing');
        return;
      }
      const vertices = drawer.graph.vertices;
      if (!Array.isArray(vertices)) {
        console.warn('annotateIsotopeCarbons: vertices is not an array');
        return;
      }

      console.log('annotateIsotopeCarbons: vertex count =', vertices.length);

      let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
      vertices.forEach(v => {
        const p = v.position;
        if (p) {
          if (p.x < minX) minX = p.x;
          if (p.x > maxX) maxX = p.x;
          if (p.y < minY) minY = p.y;
          if (p.y > maxY) maxY = p.y;
        }
      });

      const pad = 10;
      const canvas = drawer.canvasWrapper.canvas;
      const ctx = drawer.canvasWrapper.ctx;
      if (!ctx) {
        console.warn('annotateIsotopeCarbons: canvas context missing');
        return;
      }
      ctx.save();
      ctx.fillStyle = '#000';
      ctx.font = '8pt Droid Sans, sans-serif';

      vertices.forEach((v, i) => {
        const atom = v.value;
        if (!atom) return;

        const bracket = atom.bracket;
        const isotope = bracket && typeof bracket.isotope === 'number' ? bracket.isotope : null;
        if (!isotope) return;

        const p = v.position;
        if (!p) return;

        console.log('vertex', i, 'element=', atom.element, 'isotope=', isotope, 'pos=', p);

        const sx = (p.x - minX) / ((maxX - minX) || 1);
        const sy = (p.y - minY) / ((maxY - minY) || 1);
        const screenX = pad + sx * (canvas.width - 2 * pad);
        const screenY = pad + sy * (canvas.height - 2 * pad);

        ctx.fillText(String(isotope), screenX, screenY);
      });
      ctx.restore();
    }
  
    function drawMoleculeOnCanvas(smiles, canvas) {
      SmilesDrawer.parse(
        smiles,
        function (tree) {
          drawer.draw(tree, canvas.id, "light", false);
          try {
            annotateIsotopeCarbons();
          } catch (err) {
            console.error("Error annotating isotope carbons:", err);
          }
        },
        function (err) {
          console.error("SMILES parse error (" + smiles + "):", err);
        }
      );
    }
  
    function buildLinearGroup(smilesList) {
      const col = document.createElement("div");
      col.style.display = "flex";
      col.style.flexDirection = "row";
      col.style.alignItems = "center";
      col.style.gap = "8px";
  
      smilesList.forEach((smiles, i) => {
        if (!smiles || !smiles.trim()) return;
  
        if (i > 0) {
          const plus = document.createElement("div");
          plus.textContent = "+";
          plus.style.fontSize = "18px";
          plus.style.lineHeight = "1";
          col.appendChild(plus);
        }
  
        const canvas = document.createElement("canvas");
        canvas.width = drawerOptions.width;
        canvas.height = drawerOptions.height;
        canvas.id = "mol-" + Math.random().toString(36).slice(2);
        canvas.dataset.smiles = smiles;
        canvas.style.background = "#fff";
  
        col.appendChild(canvas);
      });
  
      return col;
    }
  
    function renderReactionView(reaction) {
      const substrates = (reaction.substrates || []).filter(s => s && s.trim());
      const products = (reaction.products || []).filter(s => s && s.trim());
  
      const wrapper = document.createElement("div");
      wrapper.style.display = "flex";
      wrapper.style.alignItems = "center";
      wrapper.style.justifyContent = "center";
      wrapper.style.gap = "20px";
  
      const leftCol = buildLinearGroup(substrates);
      const rightCol = buildLinearGroup(products);
  
      const arrow = document.createElement("div");
      arrow.textContent = "→";
      arrow.style.fontSize = "24px";
  
      wrapper.appendChild(leftCol);
      wrapper.appendChild(arrow);
      wrapper.appendChild(rightCol);
  
      return wrapper;
    }
  
    function renderReactionInto(container, reaction) {
      container.innerHTML = "";
      const view = renderReactionView(reaction);
      container.appendChild(view);
  
      const canvases = view.querySelectorAll("canvas[data-smiles]");
      canvases.forEach(cv => {
        const s = cv.dataset.smiles;
        if (s) {
          drawMoleculeOnCanvas(s, cv);
        }
      });
    }
  
    global.EVODEXReactionRenderer = {
      renderReactionInto
    };
  })(window);