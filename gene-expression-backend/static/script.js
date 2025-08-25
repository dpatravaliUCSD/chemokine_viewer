// Data structures
let partners = {};
let availableTissues = {};
let partnerGeneSet = new Set();

// Initialize the application
window.addEventListener('DOMContentLoaded', () => {
	initializeUI();
	bootstrapData();
});

function initializeUI() {
	const tissueSelect = document.getElementById('tissue-select');
	const gene1Select = document.getElementById('gene1-select');
	const gene2Select = document.getElementById('gene2-select');
	
	tissueSelect.addEventListener('change', handleTissueChange);
	gene1Select.addEventListener('change', handleGene1Change);
	gene2Select.addEventListener('change', handleGene2Change);
}

async function bootstrapData() {
	clearImage();
	try {
		const [tissuesResp, partnersResp] = await Promise.all([
			fetch('/api/tissues'),
			fetch('/api/interacting-partners')
		]);
		if (!tissuesResp.ok) throw new Error('Failed to load tissues');
		if (!partnersResp.ok) throw new Error('Failed to load partners');
		availableTissues = await tissuesResp.json();
		const partnersJson = await partnersResp.json();
		partners = partnersJson.interacting_partners || {};
		partnerGeneSet = new Set(Object.keys(partners));
		const availableNames = Object.entries(availableTissues)
			.filter(([name, info]) => info && info.available === true)
			.map(([name]) => name);
		populateTissueDropdown(availableNames);
	} catch (err) {
		showError(`Failed to initialize: ${err.message}`);
	}
}

function populateTissueDropdown(tissues) {
	const tissueSelect = document.getElementById('tissue-select');
	tissueSelect.innerHTML = '<option value="">Select a tissue...</option>';
	for (const tissue of tissues) {
		const option = document.createElement('option');
		option.value = tissue;
		option.textContent = tissue;
		tissueSelect.appendChild(option);
	}
}

function handleTissueChange() {
	const tissue = document.getElementById('tissue-select').value;
	const gene1Select = document.getElementById('gene1-select');
	const gene2Select = document.getElementById('gene2-select');
	
	gene1Select.innerHTML = '<option value="">Select first gene...</option>';
	gene2Select.innerHTML = '<option value="">Select second gene...</option>';
	clearImage();
	
	if (!tissue) {
		gene1Select.disabled = true;
		gene2Select.disabled = true;
		return;
	}
	
	// Populate Gene 1 with dictionary keys only
	const gene1Options = Array.from(partnerGeneSet).sort();
	for (const gene of gene1Options) {
		const opt = document.createElement('option');
		opt.value = gene;
		opt.textContent = gene;
		gene1Select.appendChild(opt);
	}
	gene1Select.disabled = gene1Options.length === 0;
	gene2Select.disabled = true;
}

function handleGene1Change() {
	const gene1 = document.getElementById('gene1-select').value;
	const gene2Select = document.getElementById('gene2-select');
	
	gene2Select.innerHTML = '<option value="">Select second gene...</option>';
	clearImage();
	
	if (!gene1) {
		gene2Select.disabled = true;
		return;
	}
	
	// Populate Gene 2 strictly from the selected key's partner list
	const candidates = (partners[gene1] || []).slice().sort();
	for (const gene of candidates) {
		const opt = document.createElement('option');
		opt.value = gene;
		opt.textContent = gene;
		gene2Select.appendChild(opt);
	}
	gene2Select.disabled = candidates.length === 0;
}

function handleGene2Change() {
	const tissue = document.getElementById('tissue-select').value;
	const gene1 = document.getElementById('gene1-select').value;
	const gene2 = document.getElementById('gene2-select').value;
	if (tissue && gene1 && gene2) {
		requestPlot(tissue, gene1, gene2);
	} else {
		clearImage();
	}
}

async function requestPlot(tissue, gene1, gene2) {
	const imageContainer = document.getElementById('image-container');
	imageContainer.innerHTML = '<div class="loading"></div>';
	try {
		const url = `/generate-plot?tissue=${encodeURIComponent(tissue)}&gene1=${encodeURIComponent(gene1)}&gene2=${encodeURIComponent(gene2)}`;
		const resp = await fetch(url);
		if (!resp.ok) {
			const text = await resp.text();
			throw new Error(text || `HTTP ${resp.status}`);
		}
		const blob = await resp.blob();
		const imgUrl = URL.createObjectURL(blob);
		const img = new Image();
		img.onload = () => {
			imageContainer.innerHTML = '';
			img.className = 'plot-image';
			img.alt = `${gene1} vs ${gene2} expression in ${tissue}`;
			imageContainer.appendChild(img);
			URL.revokeObjectURL(imgUrl);
		};
		img.onerror = () => {
			showError('Failed to render plot image');
		};
		img.src = imgUrl;
	} catch (err) {
		showError(`Plot request failed: ${err.message}`);
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