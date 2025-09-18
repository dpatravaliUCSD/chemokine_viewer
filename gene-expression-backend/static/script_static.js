// Static site script: use pre-generated images from S3/CDN
let partners = {};
let tissues = {};

window.addEventListener('DOMContentLoaded', () => {
	initUI();
	bootstrap();
});

function initUI() {
	const t = document.getElementById('tissue-select');
	const g1 = document.getElementById('gene1-select');
	const g2 = document.getElementById('gene2-select');
	t.addEventListener('change', onTissue);
	g1.addEventListener('change', onGene1);
	g2.addEventListener('change', onGene2);
}

async function bootstrap() {
	clearImage();
	try {
		const [tResp, pResp] = await Promise.all([
			fetch('/api/tissues'),
			fetch('/api/interacting-partners')
		]);
		tissues = await tResp.json();
		const pJson = await pResp.json();
		partners = pJson.interacting_partners || {};
		const names = Object.entries(tissues).filter(([n, info]) => info && info.available === true).map(([n]) => n);
		populateTissues(names);
	} catch (e) {
		showError('Failed to initialize');
	}
}

function populateTissues(list) {
	const t = document.getElementById('tissue-select');
	t.innerHTML = '<option value="">Select a tissue...</option>';
	for (const n of list) {
		const opt = document.createElement('option');
		opt.value = n; opt.textContent = n; t.appendChild(opt);
	}
}

function onTissue() {
	const g1 = document.getElementById('gene1-select');
	const g2 = document.getElementById('gene2-select');
	g1.innerHTML = '<option value="">Select first gene...</option>';
	g2.innerHTML = '<option value="">Select second gene...</option>';
	clearImage();
	const keys = Object.keys(partners).sort();
	for (const k of keys) {
		const opt = document.createElement('option');
		opt.value = k; opt.textContent = k; g1.appendChild(opt);
	}
	g1.disabled = keys.length === 0; g2.disabled = true;
}

function onGene1() {
	const g1 = document.getElementById('gene1-select').value;
	const g2 = document.getElementById('gene2-select');
	g2.innerHTML = '<option value="">Select second gene...</option>';
	clearImage();
	const candidates = (partners[g1] || []).slice().sort();
	for (const g of candidates) {
		const opt = document.createElement('option');
		opt.value = g; opt.textContent = g; g2.appendChild(opt);
	}
	g2.disabled = candidates.length === 0;
}

function onGene2() {
	const tissue = document.getElementById('tissue-select').value;
	const gene1 = document.getElementById('gene1-select').value;
	const gene2 = document.getElementById('gene2-select').value;
	if (tissue && gene1 && gene2) {
		displayImage(tissue, gene1, gene2);
	}
}

function displayImage(tissue, g1, g2) {
	const base = (window.IMAGE_BASE_URL || '').replace(/\/$/, '');
	const first = `${base}/${encodeURIComponent(tissue)}/${g1}_${g2}.png`;
	const second = `${base}/${encodeURIComponent(tissue)}/${g2}_${g1}.png`;
	const img = new Image();
	const container = document.getElementById('image-container');
	container.innerHTML = '<div class="loading"></div>';
	img.onload = () => { container.innerHTML = ''; img.className = 'plot-image'; container.appendChild(img); };
	img.onerror = () => {
		const img2 = new Image();
		img2.onload = () => { container.innerHTML = ''; img2.className = 'plot-image'; container.appendChild(img2); };
		img2.onerror = () => showError('Image not found');
		img2.src = second;
	};
	img.src = first;
}

function clearImage() {
	document.getElementById('image-container').innerHTML = '<div class="placeholder"><p>Select tissue and genes to view expression plots</p></div>';
}

function showError(msg) {
	document.getElementById('image-container').innerHTML = `<div class="error-message">${msg}</div>`;
} 