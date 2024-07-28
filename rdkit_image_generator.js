document.addEventListener("DOMContentLoaded", function() {
    const urlParams = new URLSearchParams(window.location.search);
    const evodexId = urlParams.get('id');
    if (!evodexId) {
        console.error('No EVODEX ID provided in the query string.');
        return;
    }

    fetch('../data/optimized_data.json')
        .then(response => response.json())
        .then(data => {
            renderEvodexRO(data, evodexId);
        })
        .catch(error => console.error('Error loading data:', error));
});

function renderEvodexRO(data, evodexId) {
    const evodex = data['EVODEX-RO'][evodexId];
    if (!evodex) {
        console.error('EVODEX ID not found in data:', evodexId);
        return;
    }

    document.getElementById('evodex-id').textContent = evodex.id;
    document.getElementById('evodex-title').textContent = evodex.title || 'Reaction Operator';
    document.getElementById('operator-smirks').textContent = evodex.smirks;

    // Render SMIRKS image using RDKit.js
    initRDKitModule().then(function (RDKit) {
        const reaction = RDKit.get_rxn(evodex.smirks);
        const svg = reaction.get_svg();

        const canvas = document.getElementById('reactionCanvas');
        const ctx = canvas.getContext('2d');
        const img = new Image();
        img.src = 'data:image/svg+xml;base64,' + btoa(svg);
        img.onload = function () {
            ctx.drawImage(img, 0, 0);
        };
    });

    const sourcePartialReactionsDiv = document.getElementById('source-partial-reactions');
    evodex.sources.forEach(sourceId => {
        const sourceDiv = document.createElement('div');
        const sourceLink = document.createElement('a');
        sourceLink.href = `../pages/${sourceId}.html`;
        sourceLink.textContent = sourceId;
        const sourceCanvas = document.createElement('canvas');
        sourceCanvas.width = 500;
        sourceCanvas.height = 150;
        const sourceText = document.createElement('p');
        sourceText.textContent = evodex.partial_reactions[sourceId] || 'No partial reaction data available';

        sourceDiv.appendChild(sourceLink);
        sourceDiv.appendChild(sourceCanvas);
        sourceDiv.appendChild(sourceText);
        sourcePartialReactionsDiv.appendChild(sourceDiv);

        // Render partial reaction image using RDKit.js
        initRDKitModule().then(function (RDKit) {
            const reaction = RDKit.get_rxn(evodex.partial_reactions[sourceId]);
            const svg = reaction.get_svg();

            const ctx = sourceCanvas.getContext('2d');
            const img = new Image();
            img.src = 'data:image/svg+xml;base64,' + btoa(svg);
            img.onload = function () {
                ctx.drawImage(img, 0, 0);
            };
        });
    });
}
