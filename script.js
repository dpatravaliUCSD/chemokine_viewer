// Data structure to hold tissue and gene information
const tissueData = {
    "Liver": [],
    "MC38 tumor": []
};

// Initialize the application
document.addEventListener('DOMContentLoaded', function() {
    loadTissueData();
    setupEventListeners();
});

function setupEventListeners() {
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    
    tissueSelect.addEventListener('change', handleTissueChange);
    gene1Select.addEventListener('change', handleGene1Change);
    gene2Select.addEventListener('change', handleGene2Change);
}

function loadTissueData() {
    // Define the available PNG files for each tissue
    const liverFiles = [
        "Xcr1_Xcl1.png", "Xcl1_Xcr1.png", "Tgfbr2_Tgfb3.png", "Tgfbr2_Tgfb2.png", "Tgfbr2_Tgfb1.png",
        "Tgfbr1_Tgfb3.png", "Tgfbr1_Tgfb2.png", "Tgfbr1_Tgfb1.png", "Tgfb3_Tgfbr2.png", "Tgfb3_Tgfbr1.png",
        "Tgfb2_Tgfbr2.png", "Tgfb2_Tgfbr1.png", "Tgfb1_Tgfbr2.png", "Tgfb1_Tgfbr1.png", "Il7r_Il7.png",
        "Il7_Il7r.png", "Il7_Il2rg.png", "Il6st_Il6.png", "Il6ra_Il6.png", "Il6_Il6st.png", "Il6_Il6ra.png",
        "Il4ra_Il4.png", "Il4_Il4ra.png", "Il4_Il2rg.png", "Il2rg_Il7.png", "Il2rg_Il4.png", "Il2rg_Il2.png",
        "Il2rg_Il15.png", "Il2rb_Il2.png", "Il2rb_Il15.png", "Il2ra_Il2.png", "Il2_Il2rg.png", "Il2_Il2rb.png",
        "Il2_Il2ra.png", "Il18rap_Il18.png", "Il18r1_Il18.png", "Il18_Il18rap.png", "Il18_Il18r1.png",
        "Il17rc_Il17f.png", "Il17rc_Il17a.png", "Il17ra_Il17f.png", "Il17ra_Il17a.png", "Il17f_Il17rc.png",
        "Il17f_Il17ra.png", "Il17a_Il17rc.png", "Il17a_Il17ra.png", "Il15ra_Il15.png", "Il15_Il2rg.png",
        "Il15_Il2rb.png", "Il15_Il15ra.png", "Il10rb_Il10.png", "Il10ra_Il10.png", "Il10_Il10rb.png",
        "Il10_Il10ra.png", "Ifngr2_Ifng.png", "Ifngr1_Ifng.png", "Ifng_Ifngr2.png", "Ifng_Ifngr1.png",
        "Flt3l_Flt3.png", "Flt3_Flt3l.png", "Cxcr6_Cxcl16.png", "Cxcr5_Cxcl13.png", "Cxcr4_Cxcl12.png",
        "Cxcr3_Cxcl9.png", "Cxcr3_Cxcl10.png", "Cxcl9_Cxcr3.png", "Cxcl16_Cxcr6.png", "Cxcl13_Cxcr5.png",
        "Cxcl12_Cxcr4.png", "Cxcl10_Cxcr3.png", "Cx3cr1_Cx3cl1.png", "Cx3cl1_Cx3cr1.png", "Ccr9_Ccl25.png",
        "Ccr8_Ccl8.png", "Ccr8_Ccl4.png", "Ccr8_Ccl17.png", "Ccr8_Ccl1.png", "Ccr7_Ccl19.png", "Ccr6_Ccl20.png",
        "Ccr5_Ccl8.png", "Ccr5_Ccl5.png", "Ccr5_Ccl4.png", "Ccr5_Ccl3.png", "Ccr4_Ccl5.png", "Ccr4_Ccl3.png",
        "Ccr4_Ccl22.png", "Ccr4_Ccl2.png", "Ccr4_Ccl17.png", "Ccr2_Ccl8.png", "Ccr2_Ccl7.png", "Ccr2_Ccl2.png",
        "Ccr2_Ccl12.png", "Ccr1_Ccl7.png", "Ccr1_Ccl5.png", "Ccr1_Ccl4.png", "Ccr1_Ccl3.png", "Ccl8_Ccr8.png",
        "Ccl8_Ccr5.png", "Ccl8_Ccr2.png", "Ccl7_Ccr2.png"
    ];
    
    const mc38Files = [
        "Xcr1_Xcl1.png", "Xcl1_Xcr1.png", "Tgfbr2_Tgfb3.png", "Tgfbr2_Tgfb2.png", "Tgfbr2_Tgfb1.png",
        "Tgfbr1_Tgfb3.png", "Tgfbr1_Tgfb2.png", "Tgfbr1_Tgfb1.png", "Tgfb3_Tgfbr2.png", "Tgfb3_Tgfbr1.png",
        "Tgfb2_Tgfbr2.png", "Tgfb2_Tgfbr1.png", "Tgfb1_Tgfbr2.png", "Tgfb1_Tgfbr1.png", "Il7r_Il7.png",
        "Il7_Il7r.png", "Il7_Il2rg.png", "Il6st_Il6.png", "Il6ra_Il6.png", "Il6_Il6st.png", "Il6_Il6ra.png",
        "Il4ra_Il4.png", "Il4_Il4ra.png", "Il4_Il2rg.png", "Il2rg_Il7.png", "Il2rg_Il4.png", "Il2rg_Il2.png",
        "Il2rg_Il15.png", "Il2rb_Il2.png", "Il2rb_Il15.png", "Il2ra_Il2.png", "Il2_Il2rg.png", "Il2_Il2rb.png",
        "Il2_Il2ra.png", "Il18rap_Il18.png", "Il18r1_Il18.png", "Il18_Il18rap.png", "Il18_Il18r1.png",
        "Il17rc_Il17f.png", "Il17rc_Il17a.png", "Il17ra_Il17f.png", "Il17ra_Il17a.png", "Il17f_Il17rc.png",
        "Il17f_Il17ra.png", "Il17a_Il17rc.png", "Il17a_Il17ra.png", "Il15ra_Il15.png", "Il15_Il2rg.png",
        "Il15_Il2rb.png", "Il15_Il15ra.png", "Il10rb_Il10.png", "Il10ra_Il10.png", "Il10_Il10rb.png",
        "Il10_Il10ra.png", "Ifngr2_Ifng.png", "Ifngr1_Ifng.png", "Ifng_Ifngr2.png", "Ifng_Ifngr1.png",
        "Flt3l_Flt3.png", "Flt3_Flt3l.png", "Cxcr6_Cxcl16.png", "Cxcr5_Cxcl13.png", "Cxcr4_Cxcl12.png",
        "Cxcr3_Cxcl9.png", "Cxcr3_Cxcl10.png", "Cxcl9_Cxcr3.png", "Cxcl16_Cxcr6.png", "Cxcl13_Cxcr5.png",
        "Cxcl12_Cxcr4.png", "Cxcl10_Cxcr3.png", "Cx3cr1_Cx3cl1.png", "Cx3cl1_Cx3cr1.png", "Ccr9_Ccl25.png",
        "Ccr8_Ccl8.png", "Ccr8_Ccl4.png", "Ccr8_Ccl17.png", "Ccr8_Ccl1.png", "Ccr7_Ccl19.png", "Ccr6_Ccl20.png",
        "Ccr5_Ccl8.png", "Ccr5_Ccl5.png", "Ccr5_Ccl4.png", "Ccr5_Ccl3.png", "Ccr4_Ccl5.png", "Ccr4_Ccl3.png",
        "Ccr4_Ccl22.png", "Ccr4_Ccl2.png", "Ccr4_Ccl17.png", "Ccr2_Ccl8.png", "Ccr2_Ccl7.png", "Ccr2_Ccl2.png",
        "Ccr2_Ccl12.png", "Ccr1_Ccl7.png", "Ccr1_Ccl5.png", "Ccr1_Ccl4.png", "Ccr1_Ccl3.png", "Ccl8_Ccr8.png",
        "Ccl8_Ccr5.png", "Ccl8_Ccr2.png", "Ccl7_Ccr2.png"
    ];
    
    // Process files to extract unique genes
    tissueData["Liver"] = processFiles(liverFiles);
    tissueData["MC38 tumor"] = processFiles(mc38Files);
    
    // Populate tissue dropdown
    populateTissueDropdown();
}

function processFiles(files) {
    const genes = new Set();
    const combinations = [];
    
    files.forEach(file => {
        const baseName = file.replace('.png', '');
        const parts = baseName.split('_');
        if (parts.length === 2) {
            const [gene1, gene2] = parts;
            genes.add(gene1);
            genes.add(gene2);
            combinations.push({gene1, gene2, filename: file});
        }
    });
    
    return {
        genes: Array.from(genes).sort(),
        combinations: combinations
    };
}

function populateTissueDropdown() {
    const tissueSelect = document.getElementById('tissue-select');
    
    Object.keys(tissueData).forEach(tissue => {
        const option = document.createElement('option');
        option.value = tissue;
        option.textContent = tissue;
        tissueSelect.appendChild(option);
    });
}

function handleTissueChange() {
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    const selectedTissue = tissueSelect.value;
    
    // Clear previous selections
    gene1Select.innerHTML = '<option value="">Select first gene...</option>';
    gene2Select.innerHTML = '<option value="">Select second gene...</option>';
    
    // Clear image
    clearImage();
    
    if (selectedTissue) {
        const genes = tissueData[selectedTissue].genes;
        
        // Populate only gene1 dropdown
        genes.forEach(gene => {
            const option1 = document.createElement('option');
            option1.value = gene;
            option1.textContent = gene;
            gene1Select.appendChild(option1);
        });
        
        // Enable gene1 select, keep gene2 disabled until gene1 is selected
        gene1Select.disabled = false;
        gene2Select.disabled = true;
    } else {
        // Disable gene selects
        gene1Select.disabled = true;
        gene2Select.disabled = true;
    }
}

function handleGene1Change() {
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    
    const selectedTissue = tissueSelect.value;
    const selectedGene1 = gene1Select.value;
    
    // Clear gene2 selection and image
    gene2Select.innerHTML = '<option value="">Select second gene...</option>';
    clearImage();
    
    if (selectedTissue && selectedGene1) {
        const combinations = tissueData[selectedTissue].combinations;
        
        // Find all genes that can pair with selectedGene1
        const availableGene2s = new Set();
        combinations.forEach(combo => {
            if (combo.gene1 === selectedGene1) {
                availableGene2s.add(combo.gene2);
            }
            if (combo.gene2 === selectedGene1) {
                availableGene2s.add(combo.gene1);
            }
        });
        
        // Populate gene2 dropdown with available partners
        const sortedGenes = Array.from(availableGene2s).sort();
        sortedGenes.forEach(gene => {
            const option = document.createElement('option');
            option.value = gene;
            option.textContent = gene;
            gene2Select.appendChild(option);
        });
        
        // Enable gene2 select
        gene2Select.disabled = false;
    } else {
        // Disable gene2 select
        gene2Select.disabled = true;
    }
}

function handleGene2Change() {
    const tissueSelect = document.getElementById('tissue-select');
    const gene1Select = document.getElementById('gene1-select');
    const gene2Select = document.getElementById('gene2-select');
    
    const selectedTissue = tissueSelect.value;
    const selectedGene1 = gene1Select.value;
    const selectedGene2 = gene2Select.value;
    
    if (selectedTissue && selectedGene1 && selectedGene2) {
        displayImage(selectedTissue, selectedGene1, selectedGene2);
    } else {
        clearImage();
    }
}

function displayImage(tissue, gene1, gene2) {
    const imageContainer = document.getElementById('image-container');
    
    // Show loading spinner
    imageContainer.innerHTML = '<div class="loading"></div>';
    
    // Find the matching combination
    const combinations = tissueData[tissue].combinations;
    let matchingFile = null;
    
    // Try both orders: gene1_gene2 and gene2_gene1
    const combination1 = combinations.find(combo => 
        combo.gene1 === gene1 && combo.gene2 === gene2
    );
    const combination2 = combinations.find(combo => 
        combo.gene1 === gene2 && combo.gene2 === gene1
    );
    
    matchingFile = combination1 || combination2;
    
    if (matchingFile) {
        const imagePath = `data/${tissue}/${matchingFile.filename}`;
        
        const img = new Image();
        img.onload = function() {
            imageContainer.innerHTML = '';
            img.className = 'plot-image';
            img.alt = `${gene1} vs ${gene2} expression in ${tissue}`;
            imageContainer.appendChild(img);
        };
        
        img.onerror = function() {
            showError(`Could not load image: ${imagePath}`);
        };
        
        img.src = imagePath;
    } else {
        showError(`No data available for ${gene1} and ${gene2} combination in ${tissue}`);
    }
}

function clearImage() {
    const imageContainer = document.getElementById('image-container');
    imageContainer.innerHTML = '<div class="placeholder"><p>Select tissue and genes to view expression plots</p></div>';
}

function showError(message) {
    const imageContainer = document.getElementById('image-container');
    imageContainer.innerHTML = `<div class="error-message">${message}</div>`;
} 