let thresholds = {};
let cnames = {};

$(document).ready(function() {
    google.charts.load('current', {
        callback: function () {
            loadCharts();
        },
        'packages':['corechart']}
    );
});

function loadCharts() {
    var thresh_load = $.ajax({
        url: '../results.thresh',
        dataType: 'text',
        success: function(resp) {
            var thresh_table = resp.split('\n');
            for (var i = 0; i < thresh_table.length; i++) {
                var line = thresh_table[i].split('\t')
                if (!thresholds.hasOwnProperty(line[0])) {
                    thresholds[line[0]] = {};
                }
                thresholds[line[0]][line[1]] = parseFloat(line[2]);
            }
            console.log(thresholds);
        },
        error: function(req, status, err) {
            console.log("something went wrong", status, err);
        }
    });

    var cname_load = $.ajax({
        url: '../cnames.txt',
        dataType: 'text',
        success: function(resp) {
            var rows = resp.split('\n');
            for (var i = 0; i < rows.length; i++) {
                var line = rows[i].split(" ");
                cnames[line[0]] = line[2];
            }
            console.log(cnames);
        },
        error: function(req, status, err) {
            console.log("something went wrong", status, err);
        }
    });

    const urlParams = new URLSearchParams(window.location.search);
    const consensus = urlParams.get('consensus');

    $.when(thresh_load, cname_load).then(function() {
        var matrices = Object.keys(thresholds[consensus]).sort();
        var hits = {};
        var cache_load = $.ajax({
            url: `../cache/${consensus}.json`,
            dataType: 'text',
            success: function(resp) {
                hits = JSON.parse(resp);
                console.log(hits);
            },
            error: function(req, status, err) {
                console.log("couldn't load cache");
                for (var i = 0; i < matrices.length; i++) {
                    loadHits(consensus, matrices[i]);
                }
            }
        });

        $.when(cache_load).then(function() {
            for (var i = 0; i < matrices.length; i++) {
                cacheToChart(consensus, matrices[i], hits);
            }
        });
    });
}

function cacheToChart(consensus, matrix, hits) {
    $('.charts').append(`<div class='chartWithOverlay' id='${matrix}'></div>`);
    $(`#${matrix}.chartWithOverlay`).append(`<div id='chart_${matrix}' class='chart'></div>`);
    $(`#${matrix}.chartWithOverlay`).append(`<div id='overlay_${matrix}' class='overlay'></div>`);

    var threshold = thresholds[consensus][matrix]

    var genomic_hits = new Array(10001).fill(0);
    for (var j = 0; j < hits.genomic[matrix].length; j++) {
        var hit = hits.genomic[matrix][j];
        if (hit <= 10000 && hit >= 0) {
            genomic_hits[hit]++;
        }
    }
    for (var j = genomic_hits.length - 2; j >= 0; j--) {
        genomic_hits[j] += genomic_hits[j + 1];
    }

    var benchmark_hits = new Array(10001).fill(0);
    for (var j = 0; j < hits.benchmark[matrix].length; j++) {
        var hit = hits.benchmark[matrix][j];
        if (hit <= 10000 && hit >= 0) {
            benchmark_hits[hit]++;
        }
    }
    for (var j = benchmark_hits.length - 2; j >= 0; j--) {
        benchmark_hits[j] += benchmark_hits[j + 1];
    }

    var xMax = Math.ceil(threshold + 60);
    var xMin = Math.ceil(threshold - 80);

    genomic_hits = genomic_hits.map((count, index) => [index, count]).slice(xMin, xMax);
    benchmark_hits = benchmark_hits.map((count, index) => [index, count]).slice(xMin, xMax);

    console.log(genomic_hits);
    console.log(benchmark_hits);
    drawChart(consensus, matrix, genomic_hits, benchmark_hits);
 }

function loadHits(consensus, matrix) {
    $('.charts').append(`<div class='chartWithOverlay' id='${matrix}'></div>`);
    $(`#${matrix}.chartWithOverlay`).append(`<div id='chart_${matrix}' class='chart'></div>`);
    $(`#${matrix}.chartWithOverlay`).append(`<div id='overlay_${matrix}' class='overlay'></div>`);

    var genomic_hits = []
    var benchmark_hits = []

    var tp_load = $.ajax({
        url: `../genomic_hits/${consensus}/${consensus}_${matrix}.sc`,
        dataType: 'text',
        success: function(resp) {
            genomic_hits = countHits(resp, thresholds[consensus][matrix]);
            console.log(genomic_hits);
        },
        error: function(req, status, err) {
            console.log("something went wrong", status, err);
        }
    });

    var fp_load = $.ajax({
        url: `../benchmark_hits/${consensus}/${consensus}_${matrix}.sc`,
        dataType: 'text',
        success: function(resp) {
            benchmark_hits = countHits(resp, thresholds[consensus][matrix]);
            console.log(benchmark_hits);
        },
        error: function(req, status, err) {
            console.log("something went wrong", status, err);
        }
    });

    $.when(tp_load, fp_load).then(function() {
        console.log("drawing chart");
        drawChart(consensus, matrix, genomic_hits, benchmark_hits);
    });
}

function countHits(resp, threshold) {
    const lines = resp.split("\n");
    const hits = new Array(10001).fill(0);
    var max = 0
    let regex = /^\s*(\d+)\s+\d+\.\d+\s+\d+\.\d+/
    for (var i = 0; i < lines.length; i++) {
        let match = lines[i].match(regex);
        if (match != null) {
            hit = parseInt(match[1])
            if (hit <= 10000 && hit >= 0) {
                hits[hit]++;
            }
        }
    }
    for (var i = hits.length - 2; i >= 0; i--) {
        hits[i] += hits[i + 1];
    }
    var xMax = Math.ceil(threshold + 60);
    var xMin = Math.ceil(threshold - 80);
    return hits.map((count, index) => [index, count]).slice(xMin, xMax);
}

function drawChart(consensus, matrix, genomic_hits, benchmark_hits) {
    var prehits = genomic_hits.map((elt, index) => {
        return [elt[0], elt[1],
                    benchmark_hits[index][1], elt[1] - benchmark_hits[index][1]];
    })
    var xMin = Math.ceil(thresholds[consensus][matrix] - 80);
    var xMax = Math.ceil(thresholds[consensus][matrix] + 60);

    var hits = [];
    hits.push(prehits[0])
    for (var i = 1; i < prehits.length; i++) {
        if (prehits[i][1] != prehits[i - 1][1] || prehits[i][2] != prehits[i - 1][2]) {
            hits.push(prehits[i]);
        }
    }
    console.log([['Complexity Adjusted Score', 'Cum_TP', 'Cum_FP', 'Cum_TP_Gain']].concat(hits));
    var data = google.visualization.arrayToDataTable([
        ['Complexity Adjusted Score', 'Cum_TP', 'Cum_FP', 'Cum_TP_Gain']].concat(hits));

    var options = {
        title: `${consensus} ${cnames[consensus]} Genomic(TP) vs Benchmark(FP) Cumulative scores for ${matrix} in hg38`,
        pointSize: 4,
        hAxis: {title: 'Complexity Adjusted Score',
            minValue: xMin, maxValue: xMax, direction: -1,
            viewWindow: {min: xMin, max: xMax}},
        vAxis: {title: 'Alignment Counts'},
    };

    var chart = new google.visualization.ScatterChart(document.getElementById(`chart_${matrix}`));

    chart.draw(data, options);
    drawVAxisLine(chart, thresholds[consensus][matrix], matrix);
    $(`#overlay_${matrix}`).html(`<div>Score threshold = ${thresholds[consensus][matrix]}</div>`);
}

function drawVAxisLine(chart, value, matrix) {
    var layout = chart.getChartLayoutInterface();
    var chartArea = layout.getChartAreaBoundingBox();

    var svg = chart.getContainer().getElementsByTagName('svg')[0];
    var xLoc = layout.getXLocation(value);
    svg.appendChild(createLine(xLoc, chartArea.top + chartArea.height, xLoc, chartArea.top, 'black', 4));
    document.querySelector(`#overlay_${matrix}`).style.top = chartArea.top + chartArea.height / 2;
    document.querySelector(`#overlay_${matrix}`).style.left = xLoc + 5;
}

function createLine(x1, y1, x2, y2, color, w) {
    var line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('x1', x1);
    line.setAttribute('y1', y1);
    line.setAttribute('x2', x2);
    line.setAttribute('y2', y2);
    line.setAttribute('stroke', color);
    line.setAttribute('stroke-width', w);
    return line;
}

