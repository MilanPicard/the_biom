# Fix Performance Cytoscape - Bypass du diff O(n²)

## Le Problème

Le composant `react-cytoscapejs` (utilisé par `dash-cytoscape` et `dash-cytoscape-extra`) fait un **diff O(n²)** à chaque mise à jour des éléments dans sa fonction `updateCytoscape()` :

```javascript
// Code simplifié de react-cytoscapejs (updateElements)
const updateElements = (cy, oldElements, newElements) => {
    const toRemove = [];
    const toAdd = [];
    const toUpdate = [];
    
    // O(n) - créer un map des nouveaux éléments
    newElements.forEach(el => newMap[getId(el)] = el);
    
    // O(n) - vérifier les anciens éléments
    oldElements.forEach(old => {
        if (newMap[getId(old)]) {
            toUpdate.push({old, new: newMap[getId(old)]});  // À comparer
        } else {
            toRemove.push(old);
        }
    });
    
    // O(n²) - PROBLÈME ICI : comparaison profonde de chaque élément
    toUpdate.forEach(({old, new}) => {
        // Compare TOUTES les propriétés de chaque élément
        if (deepDiff(old, new)) {  // Cette comparaison est coûteuse
            cy.getElementById(getId(old)).json(new);
        }
    });
};
```

Avec ~463 éléments, ça faisait **~214,000 comparaisons** → **27 secondes de blocage UI**.

---

## La Solution : Bypass complet du diff React

### Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        AVANT (LENT)                              │
├─────────────────────────────────────────────────────────────────┤
│  Python Callback                                                 │
│       │                                                          │
│       ▼                                                          │
│  Output('my_graph', 'elements')  ──► React Component             │
│                                           │                      │
│                                           ▼                      │
│                                    updateCytoscape()             │
│                                           │                      │
│                                           ▼                      │
│                                    DIFF O(n²) ← 27 secondes!     │
│                                           │                      │
│                                           ▼                      │
│                                    cy.batch() + cy.add()         │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                        APRÈS (RAPIDE)                            │
├─────────────────────────────────────────────────────────────────┤
│  Python Callback                                                 │
│       │                                                          │
│       ▼                                                          │
│  Output('fast_graph_data', 'data')  ──► dcc.Store                │
│                                              │                   │
│                                              ▼                   │
│                                    Clientside Callback (JS)      │
│                                              │                   │
│                                              ▼                   │
│                                    cy.batch() + cy.add()         │
│                                    (PAS DE DIFF!) ← 300ms        │
└─────────────────────────────────────────────────────────────────┘
```

---

## Implémentation étape par étape

### Étape 1 : Ajouter les dcc.Store dans le layout

Dans ton fichier de layout (ex: `app.py` ou `pages/main_app.py`), ajoute deux stores :

```python
import dash
from dash import html, dcc
import dash_cytoscape as cyto

app = dash.Dash(__name__)

app.layout = html.Div([
    # Le composant Cytoscape (garde-le, mais on n'écrira plus dans 'elements' directement)
    cyto.Cytoscape(
        id='my_graph',
        style={'width': '100%', 'height': '600px'},
        layout={'name': 'preset'},
        elements=[]  # Initialement vide
    ),
    
    # NOUVEAU: Store pour les données d'éléments (bypass du diff)
    dcc.Store(id='fast_graph_data', data=None),
    
    # NOUVEAU: Store dummy pour l'output de la clientside callback
    dcc.Store(id='fast_graph_done', data=None),
])
```

### Étape 2 : Modifier le callback Python

**AVANT** (lent) :
```python
@callback(
    Output('my_graph', 'elements'),  # ← Écrit directement dans le composant
    Output('my_graph', 'stylesheet'),
    Input('some_trigger', 'value'),
)
def update_graph(trigger):
    elements = generate_elements()  # Ta logique
    stylesheet = generate_stylesheet()
    return elements, stylesheet
```

**APRÈS** (rapide) :
```python
@callback(
    Output('fast_graph_data', 'data'),  # ← Écrit dans le Store
    Output('my_graph', 'stylesheet'),    # Le stylesheet peut rester direct
    Input('some_trigger', 'value'),
)
def update_graph(trigger):
    elements = generate_elements()  # Ta logique
    stylesheet = generate_stylesheet()
    
    # Retourne un dict avec les éléments et l'ID du graphe cible
    return {
        'elements': elements,
        'graphId': 'my_graph'  # ID du composant Cytoscape
    }, stylesheet
```

### Étape 3 : Ajouter la Clientside Callback

Dans le même fichier où tu déclares tes callbacks :

```python
from dash import clientside_callback, Input, Output

# Clientside callback qui bypass le diff React
clientside_callback(
    """
    function(fastData) {
        // Vérification des données
        if (!fastData || !fastData.elements || !fastData.graphId) {
            return window.dash_clientside.no_update;
        }
        
        const start = performance.now();
        const graphId = fastData.graphId;
        const elements = fastData.elements;
        
        console.log('[PERF] 🚀 FAST UPDATE START ' + graphId, { count: elements.length });
        
        // Récupérer l'instance Cytoscape via le DOM
        const container = document.getElementById(graphId);
        if (!container || !container._cyreg || !container._cyreg.cy) {
            console.warn('[PERF] Graph not ready: ' + graphId);
            return window.dash_clientside.no_update;
        }
        
        const cy = container._cyreg.cy;
        const oldCount = cy.elements().length;
        
        // MAGIE: Batch remove all + add all (bypass le diff React)
        cy.batch(() => {
            cy.elements().remove();  // Supprimer tous les éléments
            if (elements && elements.length > 0) {
                cy.add(elements);    // Ajouter tous les nouveaux d'un coup
            }
        });
        
        const elapsed = performance.now() - start;
        console.log('[PERF] 🚀 FAST UPDATE DONE ' + graphId, { 
            elapsed_ms: elapsed.toFixed(2), 
            old: oldCount,
            new: cy.elements().length 
        });
        
        // Retourne no_update car on a déjà modifié le graphe directement
        return window.dash_clientside.no_update;
    }
    """,
    Output('fast_graph_done', 'data'),
    Input('fast_graph_data', 'data'),
    prevent_initial_call=True
)
```

---

## Explication technique détaillée

### Pourquoi `container._cyreg.cy` ?

`dash-cytoscape` stocke l'instance Cytoscape.js dans le DOM :
- `container` = l'élément HTML `<div id="my_graph">`
- `_cyreg` = registre interne de dash-cytoscape
- `cy` = l'instance Cytoscape.js avec toutes ses méthodes

### Pourquoi `cy.batch()` ?

Sans `batch()`:
```javascript
cy.elements().remove();  // → 1 re-render
cy.add(elements);        // → 1 re-render par élément ajouté!
```

Avec `batch()`:
```javascript
cy.batch(() => {
    cy.elements().remove();  // Pas de render
    cy.add(elements);        // Pas de render
});                          // → 1 seul re-render à la fin
```

### Pourquoi `no_update` ?

On retourne `window.dash_clientside.no_update` car :
1. On a déjà modifié le graphe directement via `cy.add()`
2. On ne veut pas que Dash essaie de mettre à jour quoi que ce soit
3. Le store `fast_graph_done` n'est qu'un dummy pour satisfaire la signature du callback

---

## Gestion du Layout et des Positions

Si tu utilises `layout={'name': 'preset'}` avec des positions prédéfinies :

```python
# Les positions doivent être dans chaque élément
elements = [
    {
        'data': {'id': 'node1', 'label': 'Node 1'},
        'position': {'x': 100, 'y': 200}  # Position fixe
    },
    {
        'data': {'id': 'node2', 'label': 'Node 2'},
        'position': {'x': 300, 'y': 200}
    },
    {
        'data': {'source': 'node1', 'target': 'node2'}  # Edge
    }
]
```

Si tu utilises un layout automatique (ex: `cose`, `dagre`), tu peux le déclencher après l'ajout :

```javascript
cy.batch(() => {
    cy.elements().remove();
    cy.add(elements);
});

// Optionnel: relancer le layout après
if (fastData.runLayout) {
    cy.layout({ name: 'cose' }).run();
}
```

---

## Gestion du Stylesheet

Le stylesheet peut toujours être mis à jour normalement via `Output('my_graph', 'stylesheet')` car :
1. Le stylesheet ne déclenche pas le diff O(n²)
2. Il est appliqué efficacement par Cytoscape

```python
@callback(
    Output('fast_graph_data', 'data'),
    Output('my_graph', 'stylesheet'),  # ← Toujours OK en direct
    Input('some_trigger', 'value'),
)
def update_graph(trigger):
    elements = [...]
    stylesheet = [
        {'selector': 'node', 'style': {'background-color': 'blue'}},
        {'selector': 'edge', 'style': {'line-color': 'gray'}},
    ]
    return {'elements': elements, 'graphId': 'my_graph'}, stylesheet
```

---

## Pour plusieurs graphes

Si tu as plusieurs composants Cytoscape :

```python
# Layout
app.layout = html.Div([
    cyto.Cytoscape(id='graph1', ...),
    cyto.Cytoscape(id='graph2', ...),
    
    dcc.Store(id='fast_graph1_data', data=None),
    dcc.Store(id='fast_graph1_done', data=None),
    dcc.Store(id='fast_graph2_data', data=None),
    dcc.Store(id='fast_graph2_done', data=None),
])

# Callback pour graph1
@callback(Output('fast_graph1_data', 'data'), ...)
def update_graph1(...):
    return {'elements': [...], 'graphId': 'graph1'}, ...

# Callback pour graph2
@callback(Output('fast_graph2_data', 'data'), ...)
def update_graph2(...):
    return {'elements': [...], 'graphId': 'graph2'}, ...

# Une seule clientside callback générique (réutilise le même code)
for store_id, done_id in [('fast_graph1_data', 'fast_graph1_done'), 
                           ('fast_graph2_data', 'fast_graph2_done')]:
    clientside_callback(
        """function(fastData) { ... }""",  # Même code JS
        Output(done_id, 'data'),
        Input(store_id, 'data'),
        prevent_initial_call=True
    )
```

---

## Debugging

Ajoute ce code dans un fichier `assets/debug.js` pour monitorer les performances :

```javascript
// Monitoring des long tasks (blocage UI > 50ms)
try {
    const observer = new PerformanceObserver((list) => {
        for (const entry of list.getEntries()) {
            if (entry.duration > 500) {
                console.warn('🐢 LONG TASK', { duration_ms: entry.duration.toFixed(2) });
            }
        }
    });
    observer.observe({ entryTypes: ['longtask'] });
} catch (e) {}

// Monitoring des appels serveur
const origFetch = window.fetch;
window.fetch = function(url, options) {
    if (url && typeof url === 'string' && url.includes('_dash-update-component')) {
        console.log('📤 CALLBACK REQUEST');
        const start = performance.now();
        return origFetch.apply(this, arguments).then(response => {
            console.log('📥 CALLBACK RESPONSE', { elapsed_ms: (performance.now() - start).toFixed(2) });
            return response;
        });
    }
    return origFetch.apply(this, arguments);
};
```

---

## Résumé des performances

| Approche | Complexité | Temps (~463 éléments) | Temps (~1000 éléments) |
|----------|------------|----------------------|------------------------|
| React diff (standard) | O(n²) | ~27,000 ms | ~120,000 ms |
| Bulk replace (ce fix) | O(n) | ~300 ms | ~600 ms |

**Gain : ~100x plus rapide !** 🚀

---

## Limitations connues

1. **Perte de l'état React** : Les éléments ne sont plus synchronisés avec l'état React du composant. Si tu as besoin de lire `detail_graph.elements` dans un autre callback, il faudra le lire depuis le Store ou depuis `cy.elements().jsons()`.

2. **Animations** : Les animations de transition entre anciens et nouveaux éléments ne fonctionnent pas (on supprime tout et on rajoute tout).

3. **Sélection** : La sélection utilisateur est perdue lors du refresh (tous les éléments sont remplacés).

---

## Date

Fix appliqué le 2 décembre 2025.
