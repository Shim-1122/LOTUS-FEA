import numpy as np

from rlfea.geometry.surface import audit_body, load_body

STEP = """ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('t'),'2;1');
FILE_NAME('t','',(''),(''),'','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=COORDINATES_LIST('',4,((0.,0.,0.),(0.01,0.,0.),(0.,0.01,0.),(0.,0.,0.01)));
#2=TRIANGULATED_FACE('f',#1,4,(),$,(),((1,3,2),(1,2,4),(2,3,4),(3,1,4)));
#3=TESSELLATED_SHELL('tet',(#2),$);
#4=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));
ENDSEC;
END-ISO-10303-21;
"""


def test_ap242_tessellated_reader_units(tmp_path):
    p = tmp_path / "tet.step"
    p.write_text(STEP)
    (b,) = load_body(str(p))
    assert b.source_unit == "METRE"
    assert np.allclose(b.vertices.max(0), [10, 10, 10])      # metres -> mm
    a = audit_body(b)
    assert a["watertight"] and abs(a["signed_volume_mm3"]) == np.float64(1000 / 6) or abs(abs(a["signed_volume_mm3"]) - 1000 / 6) < 1e-9


def test_empty_step_body(tmp_path):
    p = tmp_path / "empty.step"
    p.write_text(STEP.replace("#1=COORDINATES_LIST('',4,((0.,0.,0.),(0.01,0.,0.),(0.,0.01,0.),(0.,0.,0.01)));\n#2=TRIANGULATED_FACE('f',#1,4,(),$,(),((1,3,2),(1,2,4),(2,3,4),(3,1,4)));\n#3=TESSELLATED_SHELL('tet',(#2),$);\n", ""))
    (b,) = load_body(str(p))
    assert b.empty and audit_body(b)["status"] == "EMPTY"
