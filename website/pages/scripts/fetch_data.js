// fetch_data.js

async function fetchData() {
    const response = await fetch('../data/optimized_data.json');
    const data = await response.json();
    return data;
}

function renderReactionPage(data, evodex_id) {
    const reaction = data["EVODEX-R"][evodex_id];

    if (reaction) {
        document.querySelector('#reaction-id').textContent = reaction.id;
        document.querySelector('#reaction-smirks').textContent = reaction.smirks;
        
        const sourcesContainer = document.querySelector('#reaction-sources');
        sourcesContainer.innerHTML = '';

        reaction.sources.forEach(source => {
            const sourceElement = document.createElement('div');
            sourceElement.textContent = source;
            sourcesContainer.appendChild(sourceElement);
        });

        // Render other details if necessary
        if (reaction.details) {
            document.querySelector('#reaction-natural').textContent = reaction.details.natural;
            document.querySelector('#reaction-organism').textContent = reaction.details.organism;
            document.querySelector('#reaction-protein-refs').textContent = reaction.details.protein_refs.join(', ');
            document.querySelector('#reaction-protein-db').textContent = reaction.details.protein_db;
            document.querySelector('#reaction-ec-num').textContent = reaction.details.ec_num;
        }
    } else {
        console.error(`Reaction with ID ${evodex_id} not found.`);
    }
}

document.addEventListener('DOMContentLoaded', async () => {
    const urlParams = new URLSearchParams(window.location.search);
    const evodex_id = urlParams.get('id');
    
    if (evodex_id) {
        const data = await fetchData();
        renderReactionPage(data, evodex_id);
    } else {
        console.error('No EVODEX ID provided in URL.');
    }
});
