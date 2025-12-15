// js/data_loader.js
// Lazy loader for EVODEX.2 CSV files in EVODEX/evodex/data/
//
// It does NOT load everything on startup. Instead, each "type" is loaded
// on first use and then cached in memory.
//
// Public API:
//   EVODEXData.loadType(type)    -> Promise<{ id -> entry }>
//   EVODEXData.getEntry(type,id) -> Promise<entry|null>
//   EVODEXData.listIds(type)     -> Promise<string[]>
//   EVODEXData.TYPE_CONFIG       -> config object (for debugging)

/**
 * Map logical EVODEX types to CSV files.
 *
 * Keys here are the "type" strings you will use in routes and views.
 * They must match the path segment you pass to EVODEXData (e.g. "E", "P", "Am", "Cm").
 *
 * The filenames below are expected in EVODEX/evodex/data relative to website/index.html.
 */
const TYPE_CONFIG = {
  // core reaction and partial-reaction operators
  R:  { file: "EVODEX/evodex/data/EVODEX-R.csv" },
  P:  { file: "EVODEX/evodex/data/EVODEX-P.csv" },
  E:  { file: "EVODEX/evodex/data/EVODEX-E.csv" },

  // B/C/D families and their "m" variants
  B:  { file: "EVODEX/evodex/data/EVODEX-B.csv" },
  Bm: { file: "EVODEX/evodex/data/EVODEX-Bm.csv" },
  C:  { file: "EVODEX/evodex/data/EVODEX-C.csv" },
  Cm: { file: "EVODEX/evodex/data/EVODEX-Cm.csv" },
  D:  { file: "EVODEX/evodex/data/EVODEX-D.csv" },
  Dm: { file: "EVODEX/evodex/data/EVODEX-Dm.csv" },

  // E operators and meta-operators
  Em: { file: "EVODEX/evodex/data/EVODEX-Em.csv" },

  // Formula operators
  F:  { file: "EVODEX/evodex/data/EVODEX-F.csv" },

  // Mass operators
  M:        { file: "EVODEX/evodex/data/EVODEX-M.csv" },
  M_SUBSET: { file: "EVODEX/evodex/data/EVODEX-M_mass_spec_subset.csv" },

  // D synthesis subset
  D_SYNTHESIS: { file: "EVODEX/evodex/data/EVODEX-D_synthesis_subset.csv" },

  // Additional helper data
  SELECTED: { file: "EVODEX/evodex/data/selected_reactions.csv" }
  // ubiquitous_metabolites.txt is not CSV, so leave it out here
};

/**
 * Minimal CSV parser with support for quoted fields and embedded commas.
 * Returns an array of objects keyed by header names.
 */
function parseCSV(text) {
  const lines = text.split(/\r?\n/).filter((line) => line.trim().length > 0);
  if (!lines.length) return [];

  const headers = parseCSVLine(lines[0]);
  const rows = [];

  for (let i = 1; i < lines.length; i++) {
    const line = lines[i];
    const fields = parseCSVLine(line);
    if (fields.length === 1 && fields[0] === "") continue; // skip blank

    const obj = {};
    for (let j = 0; j < headers.length; j++) {
      const key = headers[j];
      const val = fields[j] != null ? fields[j] : "";
      obj[key] = val;
    }
    rows.push(obj);
  }

  return rows;
}

/**
 * Parse a single CSV line into an array of fields.
 * Handles double quotes and embedded commas inside quotes.
 */
function parseCSVLine(line) {
  const result = [];
  let field = "";
  let inQuotes = false;

  for (let i = 0; i < line.length; i++) {
    const ch = line[i];

    if (ch === '"') {
      if (inQuotes && i + 1 < line.length && line[i + 1] === '"') {
        // Escaped quote inside quoted field
        field += '"';
        i++; // skip next quote
      } else {
        // Toggle quoting state
        inQuotes = !inQuotes;
      }
    } else if (ch === "," && !inQuotes) {
      // Field separator
      result.push(field);
      field = "";
    } else {
      field += ch;
    }
  }
  result.push(field);

  return result;
}

/**
 * Convert an array of row objects into a map keyed by "id".
 * If there is no "id" column, falls back to numeric keys.
 */
function rowsToIdMap(rows) {
  const map = {};
  if (!rows || !rows.length) return map;

  const sample = rows[0];
  const idKey = sample.hasOwnProperty("id")
    ? "id"
    : sample.hasOwnProperty("ID")
    ? "ID"
    : null;

  if (!idKey) {
    // Fallback: index by row number
    rows.forEach((row, idx) => {
      map[String(idx)] = row;
    });
    return map;
  }

  rows.forEach((row) => {
    const id = String(row[idKey]).trim();
    if (!id) return;
    map[id] = row;
  });

  return map;
}

/**
 * Internal lazy cache:
 *   cache[type]   -> { id -> entry }
 *   pending[type] -> a Promise in-flight for that type
 */
const EVODEXData = (function () {
  const cache = {};
  const pending = {};

  function normalizeType(type) {
    if (typeof type !== "string") return type;
    if (TYPE_CONFIG[type]) return type;

    const upper = type.toUpperCase();
    if (TYPE_CONFIG[upper]) return upper;

    const cap = type.charAt(0).toUpperCase() + type.slice(1);
    if (TYPE_CONFIG[cap]) return cap;

    return type;
  }

  function loadType(rawType) {
    const type = normalizeType(rawType);
    const cfg = TYPE_CONFIG[type];

    if (!cfg) {
      return Promise.reject(new Error(`EVODEXData: unknown type '${rawType}'`));
    }

    if (cache[type]) {
      return Promise.resolve(cache[type]);
    }

    if (pending[type]) {
      return pending[type];
    }

    const url = cfg.file;

    pending[type] = fetch(url)
      .then((resp) => {
        if (!resp.ok) {
          throw new Error(`EVODEXData: failed to load ${url} (status ${resp.status})`);
        }
        return resp.text();
      })
      .then((text) => {
        const rows = parseCSV(text);
        const map = rowsToIdMap(rows);
        cache[type] = map;
        return map;
      })
      .catch((err) => {
        console.error(err);
        delete pending[type];
        throw err;
      });

    return pending[type];
  }

  function getEntry(type, id) {
    return loadType(type).then((data) => {
      if (!data || typeof data !== "object") return null;
      return data[id] || null;
    });
  }

  function listIds(type) {
    return loadType(type).then((data) => {
      if (!data || typeof data !== "object") return [];
      return Object.keys(data).sort();
    });
  }

  return {
    loadType,
    getEntry,
    listIds,
    TYPE_CONFIG
  };
})();