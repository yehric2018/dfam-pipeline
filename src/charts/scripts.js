google.charts.load('current', {'packages':['corechart']});
//google.charts.setOnLoadCallback(drawChart);

let matrix = "25p35g";
let genomic_hits = [];
let benchmark_hits = [];
let tp_gain = [];
let thresholds = {};

$(document).ready(function() {
    const urlParams = new URLSearchParams(window.location.search);
    const consensus = urlParams.get('consensus');

    var thresh_load = $.ajax({
        url: '../thresholds.txt',
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

    var tp_load = $.ajax({
        url: `../genomic_hits/${consensus}/${consensus}_${matrix}.sc`,
        dataType: 'text',
        success: function(resp) {
            genomic_hits = countHits(resp);
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
            benchmark_hits = countHits(resp);
            console.log(benchmark_hits);
        },
        error: function(req, status, err) {
            console.log("something went wrong", status, err);
        }
    });

    $.when(thresh_load).then(function() {
        console.log("threshold is loaded");
    });

    $.when(thresh_load, tp_load, fp_load).then(function() {
        console.log("drawing chart");
        drawChart();
    });
});

function countHits(resp) {
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
    return hits.map((count, index) => [index, count]).slice(60, 201);
}

function drawChart(consensus, matrix, genomic_hits, benchmark_hits) {
    var prehits = genomic_hits.map((elt, index) => {
        return [elt[0], elt[1],
                    benchmark_hits[index][1], elt[1] - benchmark_hits[index][1]];
    })
    var hits = [];
    hits.push(prehits[0])
    for (var i = 1; i < prehits.length; i++) {
        if (prehits[i][1] != prehits[i - 1][1] || prehits[i][2] != prehits[i - 1][2]) {
            hits.push(prehits[i]);
        }
    }
    var data = google.visualization.arrayToDataTable([
        ['Complexity Adjusted Score', 'Cum_TP', 'Cum_FP', 'Cum_TP_Gain']].concat(hits));

    var options = {
        title: `${consensus} Genomic(TP) vs Benchmark(FP) Cumulative scores for ${matrix} in hg38`,
        pointSize: 4,
        hAxis: {title: 'Complexity Adjusted Score',
            minValue: 60, maxValue: 200, direction: -1,
            viewWindow: {min: 60, max: 200}},
        vAxis: {title: 'Alignment Counts'},
    };

    var chart = new google.visualization.ScatterChart(document.getElementById('chart_div'));

    chart.draw(data, options);
    drawVAxisLine(chart, thresholds[consensus][matrix]);
    $('.overlay').html(`<div>Score threshold = ${thresholds[consensus][matrix]}</div>`);
}

function drawVAxisLine(chart, value) {
    var layout = chart.getChartLayoutInterface();
    var chartArea = layout.getChartAreaBoundingBox();

    var svg = chart.getContainer().getElementsByTagName('svg')[0];
    var xLoc = layout.getXLocation(value);
    svg.appendChild(createLine(xLoc, chartArea.top + chartArea.height, xLoc, chartArea.top, 'black', 4));
    document.querySelector('.overlay').style.top = chartArea.top + chartArea.height / 2;
    document.querySelector('.overlay').style.left = xLoc + 5;
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

