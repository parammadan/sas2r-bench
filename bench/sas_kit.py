"""Build a kit for running the corpus in real SAS (SAS OnDemand for Academics), so every R
translation can be compared with its program's true SAS output.

    python3 bench/sas_kit.py            # writes work/sas_kit/ and work/sas_kit.zip

Kit layout (upload the whole folder to your SAS OnDemand home, e.g. ~/sas2r_kit):
    data/*.xpt               the CDISC Pilot data the programs read
    programs/<id>.sas        each corpus program, with its data URLs pointed at data/ (nothing else changed)
    setup.sas                sets &kit and the libraries the CDISC programs expect
    run_<id>.sas             runs one program and exports every dataset it created to out/<id>/*.csv
    run_all.sas              runs them all, one after another
Then download out/ and put it at work/sas_out/ here, and run:  Rscript bench/compare_sas.R <translator>
"""

import pathlib
import re
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from corpus import DATA, ROOT, fetch_all, load_corpus, program_path  # noqa: E402

KIT = ROOT / "work" / "sas_kit"

SETUP = """/* sas2r-bench kit. Set KIT to where you uploaded this folder, then run. */
%let kit = %sysget(HOME)/sas2r_kit;
options nofmterr validvarname=upcase;

/* The CDISC programs read SDTM/ADaM libraries; convert the transport files into one library. */
libname kitdata "&kit/sasdata";
%macro load(name);
  libname _x xport "&kit/data/&name..xpt";
  proc copy in=_x out=kitdata; run;
  libname _x clear;
%mend;
%load(ae) %load(adsl) %load(adae) %load(advs) %load(adpp) %load(adcm)
"""

EXPORT = """
/* Export every dataset this program left in WORK (and anything new in KITOUT) as CSV. */
%macro export_all(prog);
  %let ds=;
  proc sql noprint;
    select catx('.', libname, memname) into :ds separated by ' ' from dictionary.tables
    where libname = 'WORK' or libname = 'KITOUT';
  quit;
  %do i = 1 %to %sysfunc(countw(&ds, %str( )));
    %let d = %scan(&ds, &i, %str( ));
    proc export data=&d outfile="&kit/out/&prog/%scan(&d, 2, .).csv" dbms=csv replace; run;
  %end;
  proc datasets lib=work kill nolist; quit;
  proc datasets lib=kitout kill nolist; quit;
%mend;
"""


def localize(text):
    """Point the program's data at the kit: URL filenames to data/<name>.xpt, CDISC libraries to KITDATA."""
    text = re.sub(r'(?i)filename\s+(\w+)\s+url\s+"[^"]*/(\w+)\.xpt"\s*;',
                  lambda m: f'filename {m.group(1)} "&kit/data/{m.group(2)}.xpt";', text)
    text = re.sub(r'(?im)^\s*%include\s+"[^"]*setup\.sas"\s*;', "/* setup.sas is replaced by the kit's setup */", text)
    text = re.sub(r'(?im)^\s*libname\s+adam\s+"[^"]*"\s*;', "libname adam (kitdata);", text)
    return text


def main():
    fetch_all()
    if KIT.exists():
        shutil.rmtree(KIT)
    (KIT / "programs").mkdir(parents=True)
    (KIT / "data").mkdir()
    (KIT / "sasdata").mkdir()
    for f in sorted(DATA.glob("*.xpt")):
        shutil.copy(f, KIT / "data" / f.name)
    (KIT / "setup.sas").write_text(SETUP)
    runs = []
    for p in load_corpus():
        src = program_path(p).read_text(encoding="utf-8", errors="replace")
        (KIT / "programs" / f"{p['id']}.sas").write_text(localize(src))
        (KIT / "out" / p["id"]).mkdir(parents=True)
        run = (f'%include "%sysget(HOME)/sas2r_kit/setup.sas";\n'
               "/* adae.sas writes ADAMNEW.ADAE: give it the input library plus a separate output library */\n"
               "libname kitout \"&kit/out_lib\";\n"
               "libname sdtmnew (kitdata); libname adamnew (kitout kitdata);\n"
               f"{EXPORT}\n"
               f'%include "&kit/programs/{p["id"]}.sas";\n'
               f"%export_all({p['id']})\n")
        (KIT / f"run_{p['id']}.sas").write_text(run)
        runs.append(p["id"])
    (KIT / "out_lib").mkdir()
    (KIT / "run_all.sas").write_text(
        "/* Run every program. In SAS Studio, an error in one program does not stop the next. */\n"
        + "\n".join(f'%include "%sysget(HOME)/sas2r_kit/run_{r}.sas";' for r in runs) + "\n")
    (KIT / "README.txt").write_text(__doc__)
    shutil.make_archive(str(KIT), "zip", KIT)
    print(f"kit: {KIT} ({len(runs)} programs); zip: {KIT}.zip")


if __name__ == "__main__":
    main()
