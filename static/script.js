// Data structure to hold tissue and gene information
let tissueData = {};
let availableGenes = [];
let interactingPartners = {};
let cytokines = [];

// Initialize the application
document.addEventListener('DOMContentLoaded', function() {
    loadInitialData();
    setupEventListeners();
});

function setupEventListeners() {
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    const celltype1Select = document.getElementById('celltype1-select');
    const celltype2Select = document.getElementById('celltype2-select');
    
    tissueSelect.addEventListener('change', handleTissueChange);
    gene1Select.addEventListener('change', handleGene1Change);
    gene2Select.addEventListener('change', handleGene2Change);

    // Cell type change listeners trigger plot regeneration
    if (celltype1Select) celltype1Select.addEventListener('change', handleCellTypeChange);
    if (celltype2Select) celltype2Select.addEventListener('change', handleCellTypeChange);
}

async function loadInitialData() {
    try {
        console.log('Loading initial data (no heavy data loading)...');
        
        // Load tissue list (just names and basic info - no data loading)
        const tissueResponse = await fetch('/api/tissues');
        const tissues = await tissueResponse.json();
        
        // Load gene information (from interaction database - no data loading)
        const geneResponse = await fetch('/api/genes');
        const geneData = await geneResponse.json();
        
        // Load interacting partners data (no data loading)
        const partnersResponse = await fetch('/api/interacting-partners');
        const partnersData = await partnersResponse.json();
        
        // Store the data
        tissueData = tissues;
        availableGenes = geneData.genes;
        interactingPartners = partnersData.interacting_partners;
        cytokines = partnersData.cytokines;
        
        console.log(`Loaded ${Object.keys(tissues).length} tissue types and ${geneData.count} genes`);
        console.log(`Loaded ${cytokines.length} cytokines/chemokines with ${partnersData.total_pairs} interaction pairs`);
        console.log('Available tissues:', Object.keys(tissues));
        
        // Populate tissue dropdown
        populateTissueDropdown();
        
    } catch (error) {
        console.error('Error loading initial data:', error);
        showError('Failed to load initial data from server');
    }
}

function populateTissueDropdown() {
    const tissueSelect = document.getElementById('tissue-select');
    
    // Clear existing options except the first one
    tissueSelect.innerHTML = '<option value="">Select a tissue...</option>';
    
    Object.keys(tissueData).forEach(tissue => {
        const option = document.createElement('option');
        option.value = tissue;
        const info = tissueData[tissue];
        if (info.available) {
            option.textContent = `${tissue} (${info.file_size_mb}MB file)`;
        } else {
            option.textContent = `${tissue} (unavailable)`;
            option.disabled = true;
        }
        tissueSelect.appendChild(option);
    });
}

async function handleTissueChange() {
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    const celltype1Select = document.getElementById('celltype1-select');
    const celltype2Select = document.getElementById('celltype2-select');
    const selectedTissue = tissueSelect.value;
    
    // Clear previous selections
    gene1Select.innerHTML = '<option value="">Select cytokine/chemokine...</option>';
    gene2Select.innerHTML = '<option value="">Select receptor/partner...</option>';

    // Reset celltype dropdowns
    if (celltype1Select) {
        celltype1Select.innerHTML = '<option value="">Select cell type (optional)...</option>';
        celltype1Select.disabled = true;
    }
    if (celltype2Select) {
        celltype2Select.innerHTML = '<option value="">Select cell type (optional)...</option>';
        celltype2Select.disabled = true;
    }
    
    // Clear image
    clearImage();
    
    if (selectedTissue) {
        console.log(`Loading metadata for ${selectedTissue}...`);
        
        try {
            // NOW load tissue metadata (only when tissue is selected)
            const metadataResponse = await fetch(`/api/tissue-metadata?tissue=${encodeURIComponent(selectedTissue)}`);
            const metadata = await metadataResponse.json();
            
            if (!metadataResponse.ok) {
                console.error('Error loading tissue metadata:', metadata.error);
                showError(`Error loading ${selectedTissue}: ${metadata.error}`);
                return;
            }
            
            console.log(`${selectedTissue} metadata loaded:`, metadata);
            
            // Populate gene1 dropdown with cytokines/chemokines
            cytokines.forEach(cytokine => {
                const option1 = document.createElement('option');
                option1.value = cytokine;
                option1.textContent = cytokine;
                gene1Select.appendChild(option1);
            });
            
            // Enable gene1 select, keep gene2 disabled until gene1 is selected
            gene1Select.disabled = false;
            gene2Select.disabled = true;
            
        } catch (error) {
            console.error('Error loading tissue metadata:', error);
            showError(`Failed to load ${selectedTissue} metadata`);
        }
    } else {
        // Disable gene selects
        gene1Select.disabled = true;
        gene2Select.disabled = true;
    }
}

async function handleGene1Change() {
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    const celltype1Select = document.getElementById('celltype1-select');
    const celltype2Select = document.getElementById('celltype2-select');
    
    const selectedTissue = tissueSelect.value;
    const selectedGene1 = gene1Select.value;
    
    // Clear gene2 selection and image
    gene2Select.innerHTML = '<option value="">Select receptor/partner...</option>';
    clearImage();
    
    // Reset celltype2 until gene2 is chosen
    if (celltype2Select) {
        celltype2Select.innerHTML = '<option value="">Select cell type (optional)...</option>';
        celltype2Select.disabled = true;
    }

    // Load cell types for gene1 if both tissue and gene1 are selected
    if (celltype1Select && selectedTissue && selectedGene1) {
        await populateCellTypeDropdown(celltype1Select, selectedTissue, selectedGene1);
    }
    
    if (selectedTissue && selectedGene1 && interactingPartners[selectedGene1]) {
        // Populate gene2 dropdown with only the receptors/partners for the selected cytokine
        const partners = interactingPartners[selectedGene1];
        partners.forEach(partner => {
            const option = document.createElement('option');
            option.value = partner;
            option.textContent = partner;
            gene2Select.appendChild(option);
        });
        
        console.log(`Selected ${selectedGene1}, available partners:`, partners);
        
        // Enable gene2 select
        gene2Select.disabled = false;
    } else {
        // Disable gene2 select
        gene2Select.disabled = true;
    }
}

async function handleGene2Change() {
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    const celltype2Select = document.getElementById('celltype2-select');
    
    const selectedTissue = tissueSelect.value;
    const selectedGene1 = gene1Select.value;
    const selectedGene2 = gene2Select.value;
    
    // Load cell types for gene2 if tissue and gene2 are selected
    if (celltype2Select && selectedTissue && selectedGene2) {
        await populateCellTypeDropdown(celltype2Select, selectedTissue, selectedGene2);
    }

    if (selectedTissue && selectedGene1 && selectedGene2) {
        console.log(`Generating plot for ${selectedGene1} -> ${selectedGene2} in ${selectedTissue}`);
        displayImage(selectedTissue, selectedGene1, selectedGene2);
    } else {
        clearImage();
    }
}

async function populateCellTypeDropdown(selectElement, tissue, gene) {
    try {
        console.log(`Loading cell types for ${gene} in ${tissue}...`);
        
        // Show loading state
        selectElement.innerHTML = '<option value="">Loading cell types...</option>';
        selectElement.disabled = true;
        
        const response = await fetch(`/api/cell-types?tissue=${encodeURIComponent(tissue)}&gene=${encodeURIComponent(gene)}`);
        const data = await response.json();
        
        if (!response.ok) {
            console.error('Error loading cell types:', data.error);
            selectElement.innerHTML = '<option value="">Error loading cell types</option>';
            return;
        }
        
        // Clear and populate with cell types
        selectElement.innerHTML = '<option value="">Select cell type (optional)...</option>';
        
        data.cell_types.forEach(cellType => {
            const option = document.createElement('option');
            option.value = cellType;
            option.textContent = cellType;
            selectElement.appendChild(option);
        });
        
        // Enable the dropdown
        selectElement.disabled = false;
        
        console.log(`Loaded ${data.cell_types.length} cell types for ${gene}: ${data.cell_types.join(', ')}`);
        console.log(`${data.total_expressing_cells}/${data.total_cells} cells express ${gene}`);
        
    } catch (error) {
        console.error('Error loading cell types:', error);
        selectElement.innerHTML = '<option value="">Error loading cell types</option>';
        selectElement.disabled = true;
    }
}

function handleCellTypeChange() {
    // Regenerate plot when cell type selection changes
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    
    const selectedTissue = tissueSelect.value;
    const selectedGene1 = gene1Select.value;
    const selectedGene2 = gene2Select.value;

    if (selectedTissue && selectedGene1 && selectedGene2) {
        console.log('Cell type filter changed, regenerating plot...');
        displayImage(selectedTissue, selectedGene1, selectedGene2);
    }
}

function displayImage(tissue, gene1, gene2) {
    const imageContainer = document.getElementById('image-container');
    const celltype1Select = document.getElementById('celltype1-select');
    const celltype2Select = document.getElementById('celltype2-select');
    
    // Get selected cell types
    const celltype1 = celltype1Select ? celltype1Select.value : '';
    const celltype2 = celltype2Select ? celltype2Select.value : '';
    
    // Show loading spinner
    imageContainer.innerHTML = '<div class="loading"></div>';
    
    // Build the API URL for dynamic plot generation with cell type parameters
    let apiUrl = `/generate-plot?tissue=${encodeURIComponent(tissue)}&gene1=${encodeURIComponent(gene1)}&gene2=${encodeURIComponent(gene2)}`;
    
    if (celltype1) {
        apiUrl += `&celltype1=${encodeURIComponent(celltype1)}`;
    }
    if (celltype2) {
        apiUrl += `&celltype2=${encodeURIComponent(celltype2)}`;
    }
    
    console.log('Plot URL:', apiUrl);
    
    // Create image element and load from API
    const img = new Image();
    img.onload = function() {
        imageContainer.innerHTML = '';
        img.className = 'plot-image';
        let altText = `${gene1} vs ${gene2} expression in ${tissue}`;
        if (celltype1 || celltype2) {
            const filters = [];
            if (celltype1) filters.push(`${gene1}: ${celltype1}`);
            if (celltype2) filters.push(`${gene2}: ${celltype2}`);
            altText += ` (filtered by ${filters.join(', ')})`;
        }
        img.alt = altText;
        imageContainer.appendChild(img);
    };
    
    img.onerror = function() {
        showError(`Could not generate plot for ${gene1} and ${gene2} in ${tissue}. Please check if the genes exist in the dataset.`);
    };
    
    // Load the dynamically generated plot
    img.src = apiUrl;
}

function clearImage() {
    const imageContainer = document.getElementById('image-container');
    imageContainer.innerHTML = '<div class="placeholder"><p>Select tissue and genes to view expression plots</p></div>';
}

function showError(message) {
    const imageContainer = document.getElementById('image-container');
    imageContainer.innerHTML = `<div class="error-message">${message}</div>`;
} 