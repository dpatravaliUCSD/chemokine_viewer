# Spatial Gene Expression Viewer

A web-based tool for visualizing spatial gene expression patterns across different tissues.

## Features

- Interactive tissue selection (SI, Liver, MC38 tumor)
- Gene pair selection with dynamic dropdowns
- Real-time image display of spatial expression plots
- Responsive design for desktop and mobile
- GitHub Pages compatible
**Access site** at: `https://amonell.github.io/chemokine_viewer`

## Data Structure

```
data/
├── Liver/
│   ├── Gene1_Gene2.png
│   ├── Gene2_Gene1.png
│   └── ...
└── MC38 tumor/
    ├── Gene1_Gene2.png
    ├── Gene2_Gene1.png
    └── ...
```

## Supported Gene Types

- Chemokines and receptors (Ccl/Ccr, Cxcl/Cxcr)
- Interleukins and receptors (Il/Ilr)
- Growth factors (Tgfb/Tgfbr)
- Other cytokines and signaling molecules

## Browser Compatibility

- Chrome (recommended)
- Firefox
- Safari
- Edge

## Technical Details

- Pure HTML/CSS/JavaScript (no frameworks required)
- Client-side only (perfect for GitHub Pages)
- Responsive grid layout
- Image lazy loading with error handling 
