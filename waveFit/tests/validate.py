#!/usr/bin/env python3
"""End-to-end checks against independent ROOT objects and known pulse truth."""
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import tempfile

import ROOT

ROOT.gROOT.SetBatch(True)
BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'output' / 'tests'
OUT.mkdir(parents=True, exist_ok=True)
RUN = Path(tempfile.mkdtemp(prefix='validation_', dir=OUT))
EXE = BASE / 'wave_reader'
results = []


def run(file, *args, good=True):
    command = [str(EXE), str(file), *map(str, args)]
    p = subprocess.run(command, cwd=BASE, text=True, capture_output=True)
    (RUN / f'command_{len(list(RUN.glob("command_*.log"))):03d}.log').write_text(
        ' '.join(command) + '\n' + p.stdout + p.stderr)
    assert (p.returncode == 0) == good, p.stdout[-3000:] + p.stderr[-3000:]
    return p


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def graph_points(g):
    return [(g.GetPointX(i), g.GetPointY(i)) for i in range(g.GetN())]


def raw_waves(pad):
    found = []
    for o in pad.GetListOfPrimitives():
        if o.InheritsFrom('TPad'):
            found.extend(raw_waves(o))
        elif o.InheritsFrom('TGraph'):
            found.append((o.GetName(), graph_points(o)))
        elif o.InheritsFrom('TH1') and o.GetDimension() == 1:
            if o.GetName() == 'hframe' and o.GetEntries() == 0 and o.GetSumOfWeights() == 0:
                continue
            found.append((o.GetName(), [(o.GetBinCenter(i), o.GetBinContent(i)) for i in range(1, o.GetNbinsX()+1)]))
    return found


# Provided data: all extracted samples, metadata, output graphs and entry counts.
for name, ncanvas in [('wave_242Fm', 22), ('wave_246Fm', 8)]:
    source = BASE / 'test' / f'{name}.root'
    before = digest(source)
    out = RUN / name
    run(source, '--inspect', '--list', '--analyze-all', '--export-waveforms', '--output-dir', out)
    original = ROOT.TFile.Open(str(source))
    analysis = ROOT.TFile.Open(str(out / f'{name}_analysis.root'))
    tree = analysis.Get('WaveAnalysis')
    expected = {}
    for key in original.GetListOfKeys():
        c = key.ReadObj()
        for index, (obj, points) in enumerate(raw_waves(c)):
            expected[(c.GetName(), index)] = (obj, points)
    assert original.GetNkeys() == ncanvas
    assert tree.GetEntries() == len(expected)
    assert len(list((out / 'waveforms').glob('*.csv'))) == len(expected)
    count = 0
    for row in tree:
        obj, points = expected[(str(row.canvasName), row.objectIndex)]
        assert str(row.objectName) == obj
        assert list(zip(row.x, row.y)) == points
        assert row.chainNum >= 1 and row.run >= 0 and row.mapEntry >= 0
        assert str(row.fitStatus) == 'not_requested'
        assert math.isnan(row.deltaT)
        ident = f'{row.canvasName}_{obj}_w{row.objectIndex}'
        csv_points = [(float(r['x']), float(r['y'])) for r in csv.DictReader((out/'waveforms'/f'{ident}.csv').open())]
        assert csv_points == points
        assert graph_points(analysis.Get(f'{ident}/raw_waveform')) == points
        count += len(points)
    csv_rows = list(csv.DictReader((out/f'{name}_analysis.csv').open()))
    assert len(csv_rows) == len(expected)
    analysis.Close()
    original.Close()
    assert before == digest(source)
    results.append({'test': name, 'canvases': ncanvas, 'waveforms': len(expected), 'exact_points': count, 'input_unchanged': True})


# Synthetic nested pads, inheritance, duplicate graph names, TH1, multiple key
# cycles, missing metadata, TH2 exclusion, empty pads, and non-data primitives.
fixture = RUN / 'structure.root'
f = ROOT.TFile(str(fixture), 'RECREATE')
c = ROOT.TCanvas('ChainNum_41', 'ChainNum_41: run = 7, map_entry = 1234567890123')
c.Divide(1, 2)
p = c.cd(1)
inner = ROOT.TPad('inner', 'nested pad', 0, 0, 1, 1)
inner.Draw(); inner.cd()
g = ROOT.TGraphErrors(16)
g.SetName('same')
for i in range(16): g.SetPoint(i, i*0.25, i*i-2)
g.Draw('AL')
line = ROOT.TLine(0, 1, 1, 2); line.Draw()
p = c.cd(2)
a = ROOT.TGraphAsymmErrors(16); a.SetName('same')
for i in range(16): a.SetPoint(i, i*0.25, 50-i)
a.Draw('AL')
c.Write(); c.Write()  # newest cycle only
sub = f.mkdir('nested'); sub.cd()
c2 = ROOT.TCanvas('unlabelled', 'no event metadata')
h = ROOT.TH1D('hist', 'Histogram waveform', 16, -2, 2); h.SetDirectory(0)
for i in range(1, 17): h.SetBinContent(i, i+0.5)
h.Draw(); c2.Write()
c3 = ROOT.TCanvas('empty', 'empty pad'); c3.Write()
c4 = ROOT.TCanvas('two_dimensional', 'not a waveform')
h2 = ROOT.TH2D('hist2', '2D histogram', 4, 0, 4, 4, 0, 4); h2.SetDirectory(0)
h2.Draw(); c4.Write()
f.Close()
struct_out = RUN/'structure_output'
p = run(fixture, '--inspect', '--export-waveforms', '--analyze-all', '--draw-all', '--output-dir', struct_out)
assert '[TGraphErrors]' in p.stdout and '[TGraphAsymmErrors]' in p.stdout and '[TLine]' in p.stdout
assert 'Canvases=4 selected_waveforms=3' in p.stdout
f = ROOT.TFile.Open(str(struct_out/'structure_analysis.root'))
for row in f.Get('WaveAnalysis'):
    if row.chainNum == 41:
        assert row.run == 7 and row.mapEntry == 1234567890123
    else:
        assert row.chainNum == -1 and row.run == -1 and row.mapEntry == -1
f.Close()
assert len(list((struct_out/'drawings').glob('*.png'))) == 3
assert len(list((struct_out/'waveforms').glob('*.csv'))) == 3
points = [(float(r['x']), float(r['y'])) for r in csv.DictReader((struct_out/'waveforms'/'nested_unlabelled_hist_w0.csv').open())]
assert points == [(h.GetBinCenter(i), i+0.5) for i in range(1, 17)]
run(fixture, '--export-waveforms', '--output-dir', struct_out, good=False)
run(fixture, '--analyze', 999, good=False)
run(fixture, '--smooth', 2, good=False)
run(fixture, '--object-pattern', '[', good=False)
run(fixture, '--fit', 41, good=False)
run(fixture, '--export-waveforms', '--object-pattern', 'no_such_wave', good=False)
run(RUN/'missing.root', '--inspect', good=False)
results.append({'test': 'structure_and_cli', 'status': 'passed'})


# Analytic Gaussian pulse: rising derivative maximum is at template t=0.
# Build truth by independent linear interpolation of known template samples.
template = RUN/'truth_template.root'
f = ROOT.TFile(str(template), 'RECREATE')
tg = ROOT.TGraph(401); tg.SetName('pulse_template')
ty = []
for i in range(401):
    x = i-80
    y = math.exp(-0.5*((x-10)/10)**2)
    ty.append(y); tg.SetPoint(i, x, y)
tg.Write(); ROOT.TNamed('x_unit', 'original_x').Write(); f.Close()


def pulse(x):
    j = x+80
    if j < 0: return 0
    if j >= 400: return ty[-1]
    i = int(j)
    return ty[i] + (ty[i+1]-ty[i])*(j-i)


source = RUN/'pulses.root'
f = ROOT.TFile(str(source), 'RECREATE')
keep = []
rng = random.Random(913)
for chain in range(1, 6):
    c = ROOT.TCanvas(f'ChainNum_{chain}', f'ChainNum_{chain}: run=2, map_entry={100+chain}')
    g = ROOT.TGraph(400); g.SetName('detector')
    for i in range(400):
        y = 1000*pulse(i-80.25)
        if chain in (2, 3): y += 650*pulse(i-80.25-35.4)
        if chain == 3: y = -y
        if chain == 4: y = 0
        if chain == 5: y = 1000*pulse(i-80.25) + 650*pulse(i-80.25-0.4)
        g.SetPoint(i, i, 8192+y+rng.gauss(0, 0.5))
    g.Draw('AL'); c.Write(); keep.extend([c, g])
f.Close()
fitout = RUN/'fits'
run(source, '--analyze-all', '--template', template, '--t1', 80.25, '--max-delta', 80, '--output-dir', fitout)
rows = list(csv.DictReader((fitout/'pulses_analysis.csv').open()))
assert len(rows) == 5
for i in (1, 2):
    row = rows[i]
    assert abs(float(row['deltaT'])-35.4) < 0.03, row
    assert abs(float(row['A1'])-1000) < 2, row
    assert abs(float(row['A2'])-650) < 2, row
    assert float(row['deltaBIC']) > 100, row
    assert row['fitStatus'] == 'pileup_candidate', row
assert int(rows[2]['polarity']) == -1
assert float(rows[0]['deltaBIC']) < 10, rows[0]
assert int(rows[3]['nPulse']) == 0, rows[3]
# Delta=0 must remain finite despite coincident pulse columns.
f = ROOT.TFile.Open(str(fitout/'pulses_analysis.root'))
scan = f.Get('ChainNum_2_detector_w0/chi2_vs_deltaT')
assert scan.GetPointX(0) == 0 and math.isfinite(scan.GetPointY(0))
power = f.Get('ChainNum_2_detector_w0/matched_filter_power')
assert power.GetN() == scan.GetN()
for i in range(scan.GetN()):
    assert abs(power.GetPointY(i)-max(0, float(rows[1]['chi2_1'])-scan.GetPointY(i))) < 1e-6
for row in f.Get('WaveAnalysis'):
    ident = f'ChainNum_{row.chainNum}_detector_w0'
    for name in ('best_single_fit', 'best_double_fit', 'residual_single', 'residual_double'):
        assert f.Get(f'{ident}/{name}').GetN() == 400
f.Close()
results.append({'test': 'known_double_pulses', 'deltaT_truth': 35.4,
                'positive_deltaT': float(rows[1]['deltaT']), 'negative_deltaT': float(rows[2]['deltaT']),
                'positive_A1': float(rows[1]['A1']), 'positive_A2': float(rows[1]['A2'])})

# Build a measured template, run automatic timing and validate unit conversion.
tout = RUN/'built_template'
run(source, '--make-template', 1, '--output-dir', tout)
f = ROOT.TFile.Open(str(tout/'pulse_template.root'))
assert f.Get('pulse_template').GetN() > 100
assert abs(max(y for _, y in graph_points(f.Get('pulse_template')))-1) < 1e-12
f.Close()
autoout = RUN/'auto_fits'
run(source, '--analyze', 2, '--template', tout/'pulse_template.root', '--max-delta', 80, '--output-dir', autoout)
auto = next(csv.DictReader((autoout/'pulses_analysis.csv').open()))
assert abs(float(auto['deltaT'])-35.4) < 1.5, auto
orderout = RUN/'option_order'
run(source, '--t1', 80.25, '--fit', 2, '--template', template, '--max-delta', 0, '--output-dir', orderout)
order = next(csv.DictReader((orderout/'pulses_analysis.csv').open()))
assert float(order['t1']) == 80.25 and float(order['deltaT']) == 0
nsout = RUN/'ns'
run(source, '--analyze', 2, '--sample-period-ns', 10, '--output-dir', nsout)
ns = next(csv.DictReader((nsout/'pulses_analysis.csv').open()))
assert ns['xUnit'] == 'ns'
run(source, '--analyze', 2, '--sample-period-ns', 10, '--template', template, '--output-dir', RUN/'bad_units', good=False)
results.append({'test': 'template_and_units', 'auto_deltaT': float(auto['deltaT']), 'status': 'passed'})

report = {'status': 'passed', 'results': results, 'output': str(RUN)}
(RUN/'validation.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
