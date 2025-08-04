# Spatial Gene Expression Viewer

A web-based tool for visualizing spatial gene expression patterns across different tissues.

## Features

- Interactive tissue selection (SI, Liver, MC38 tumor)
- Gene pair selection with dynamic dropdowns
- Real-time image display of spatial expression plots
- Responsive design for desktop and mobile
- GitHub Pages compatible

## Quick Deployment to GitHub Pages

1. **Create a new repository** on GitHub
2. **Upload these files** to your repository:
   - `index.html`
   - `styles.css`
   - `script.js`
   - `data/` folder (with all tissue subfolders and PNG files)

3. **Enable GitHub Pages**:
   - Go to repository Settings
   - Scroll to Pages section
   - Select "Deploy from a branch"
   - Choose "main" branch and "/ (root)" folder
   - Click Save

4. **Access your site** at: `https://yourusername.github.io/repositoryname`

## Usage

1. Select a tissue from the dropdown
2. Choose two genes from the available options
3. View the spatial expression plot
4. The system automatically finds matching gene combinations regardless of order

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