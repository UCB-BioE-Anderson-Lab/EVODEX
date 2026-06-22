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
    <section>
      <h2>About EVODEX.2</h2>
      <p>
        EVODEX is a framework for understanding how enzymes transform molecules.
        It helps researchers figure out which products can be formed from a given starting molecule,
        determining the mechanistic plausibility of potential transformations,
        and interpreting mass spectrometry data.
        This is useful in areas like metabolic engineering, drug discovery,
        interpreting lab data, and enabling AI to reason about biochemistry.
      </p>
      <p>
        EVODEX.2 is a release of the EVODEX (Enzymatic Validation Operator Differentiation and Extraction) framework that abstracts enzymatic reactions derived from the BRENDA database
        to generate a family of named operators.
        These include representations of full and partial reactions, mass and formula changes,
        and mechanistic patterns derived from σ and π bonding interactions.
        Additionally, it provides algorithms for using the operators.
        For information on using EVODEX, visit
        <a href="https://github.com/UCB-BioE-Anderson-Lab/EVODEX" target="_blank" rel="noopener noreferrer">our GitHub repository</a>.
      </p>
      <p>
        The “.2” indicates that this is version 2 of EVODEX. All operator identifiers are namespaced with this version number (e.g., EVODEX.2-R12), so they can be uniquely referenced and distinguished from past and future versions.
      </p>
      <p>
        For more technical information on EVODEX, see our publication:
        Conceicao, L. L.; Lin, H.; Tai, Y.-C.; Xue, G.; Zhang, H.; Zhang, C.; Anderson, J. C.
        EVODEX: A Cheminformatics Framework for Mechanistic Abstraction and Prediction of Enzymatic Reactions.
        <i>J. Chem. Inf. Model.</i> 2026, published online April 2.
        <a href="https://doi.org/10.1021/acs.jcim.5c03163" target="_blank" rel="noopener noreferrer">https://doi.org/10.1021/acs.jcim.5c03163</a>.
        PMID: 41925163.
      </p>
    </section>

    <section>
      <h2>Google Colab Notebooks</h2>
      <p>
        These interactive notebooks demonstrate how EVODEX.2 operators can be applied to predict biosynthetic pathways,
        evaluate proposed reactions, and interpret mass spectrometry data. Each demo walks through real examples using
        operators mined from the EVODEX.2 dataset.
      </p>
      <ul>
        <li><a href="https://colab.research.google.com/drive/16liT8RhMCcRzXa_BVdYX7xgbgVAWK4tA" target="_blank" rel="noopener noreferrer">EVODEX Synthesis Demo</a></li>
        <li><a href="https://colab.research.google.com/drive/1IvoaXjtnu7ZSvot_1Ovq3g-h5IVCdSn4" target="_blank" rel="noopener noreferrer">EVODEX Evaluation Demo</a></li>
        <li><a href="https://colab.research.google.com/drive/1CV5HM9lBy-U-J6nLqBlO6Y1WtCFWP8rX" target="_blank" rel="noopener noreferrer">EVODEX Mass Spec Demo</a></li>
      </ul>
    </section>

    <section>
      <h2>EVODEX.2 Operator Glossary</h2>

      <table class="evodex-glossary">
        <thead>
          <tr>
            <th>Operator</th>
            <th>Description</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><a href="#/type/R">R</a></td>
            <td>Full atom-mapped reaction</td>
          </tr>
          <tr>
            <td><a href="#/type/P">P</a></td>
            <td>Partial reaction based on shared atom maps</td>
          </tr>
          <tr>
            <td><a href="#/type/F">F</a></td>
            <td>Formula difference between reactants/products</td>
          </tr>
          <tr>
            <td><a href="#/type/M">M</a></td>
            <td>Mass difference for MS analysis</td>
          </tr>

          <tr>
            <td>
              <div>A</div>
              <div>Am</div>
            </td>
            <td>
              <div>Abstract. Reaction centers devoid of atom identity.</div>
              <div>Matched.</div>
              <div><em>Not computed in EVODEX.2, but supported by the operator extraction code.</em></div>
            </td>
          </tr>

          <tr>
            <td>
              <div><a href="#/type/B">B</a></div>
              <div><a href="#/type/Bm">Bm</a></div>
            </td>
            <td>
              <div>Bond-Changing. Includes reaction centers (connectivity changes or stereochemical changes).</div>
              <div>Matched.</div>
            </td>
          </tr>

          <tr>
            <td>
              <div><a href="#/type/C">C</a></div>
              <div><a href="#/type/Cm">Cm</a></div>
            </td>
            <td>
              <div>Covalent Shell. Includes first σ shell around bond changing atoms.</div>
              <div>Matched.</div>
            </td>
          </tr>

          <tr>
            <td>
              <div><a href="#/type/D">D</a></div>
              <div><a href="#/type/Dm">Dm</a></div>
            </td>
            <td>
              <div>Delocalized Domain. Includes π shell around the σ shell.</div>
              <div>Matched.</div>
            </td>
          </tr>

          <tr>
            <td>
              <div><a href="#/type/E">E</a></div>
              <div><a href="#/type/Em">Em</a></div>
            </td>
            <td>
              <div>Electronic Domain. Includes second σ shell around the π shell.</div>
              <div>Matched.</div>
            </td>
          </tr>

          <tr>
            <td><a href="#/type/D_SYNTHESIS">D_SYNTHESIS</a></td>
            <td>D synthesis subset</td>
          </tr>
          <tr>
            <td><a href="#/type/M_SUBSET">M_SUBSET</a></td>
            <td>M mass spec subset</td>
          </tr>
        </tbody>
      </table>
    </section>

    <section>
      <h2>Focused Subsets for Synthesis and Mass Spectrometry</h2>
      <p>
        These focused subsets of the EVODEX.2 dataset support specific EVODEX use cases.
        Both subsets are constructed by restricting to reactions that can be represented as a single substrate to a single product.
        The D synthesis subset presents EVODEX-D operators for this setting. In principle, the same restriction could be applied to
        any abstraction family, but we present D because it balances generality with precise capture of reactivity patterns, while
        E operators are often over-specified for this application. The mass spectrometry subset is derived from EVODEX-M and
        similarly restricts to single-substrate to single-product transformations, yielding mass differences that plausibly relate
        two LCMS peaks through a substrate/product relationship.
      </p>
      <ul>
        <li><a href="#/type/D_SYNTHESIS">Synthesis Subset of EVODEX-D</a></li>
        <li><a href="#/type/M_SUBSET">Mass Spectrometry Subset of EVODEX-M</a></li>
      </ul>
    </section>

    <section>
      <h2>Credits and Contact</h2>
      <p>
        Developed by the Anderson Lab, Department of Bioengineering, UC Berkeley.<br>
        For more information, contact: <a href="mailto:jcanderson@berkeley.edu">jcanderson@berkeley.edu</a>
      </p>
      <p>
        Please cite our work when using EVODEX:<br>
        Conceicao, L. L.; Lin, H.; Tai, Y.-C.; Xue, G.; Zhang, H.; Zhang, C.; Anderson, J. C.
        EVODEX: A Cheminformatics Framework for Mechanistic Abstraction and Prediction of Enzymatic Reactions.
        <i>J. Chem. Inf. Model.</i> 2026, published online April 2.
        <a href="https://doi.org/10.1021/acs.jcim.5c03163" target="_blank" rel="noopener noreferrer">https://doi.org/10.1021/acs.jcim.5c03163</a>.
        PMID: 41925163.
      </p>
      <p>EVODEX was derived from the following resources:</p>
      <ul>
        <li><a href="https://www.brenda-enzymes.org/" target="_blank" rel="noopener noreferrer">BRENDA: The Comprehensive Enzyme Information System</a></li>
        <li>
          EnzymeMap: Heid E, Probst D, Green WH, Madsen GKH.
          <i>EnzymeMap: curation, validation and data-driven prediction of enzymatic reactions.</i>
          <a href="https://doi.org/10.1039/d3sc02048g" target="_blank" rel="noopener noreferrer">Chem Sci. 2023;14(48):14229–14242</a>.
          PMID: 38098707; PMCID: PMC10718068.
        </li>
      </ul>
    </section>
  `;

  app.innerHTML = card("Home", introHtml);
}