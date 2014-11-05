import re
import math
import sys
import time
import os
import glob
import os.path
import logging
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import scipy.stats
import h5py
import matplotlib.pyplot as plt
import quickpen
import penplot
import time
from default_parser import DefaultArgumentParser

logger=logging.getLogger(__file__)

def metadata(h5f, h5i):
    info={}
    attrs=h5f['/trajectory'].attrs
    info['CompileTime']=attrs['COMPILETIME'].tostring().decode("utf-8")
    info['Makefile']=attrs['CONFIG'].tostring().decode("utf-8")
    gitrepo=attrs['VERSION'].tostring().decode("utf-8")
    elapsed_ns=0
    try:
        print(list(attrs.keys()))
        #print(attrs['totalelapsedwall'][0].decode()["utf-8"])
        elapsed_ns=float(attrs['totalelapsedwall'][0])
        print(elapsed_ns)
    except:
        pass
    info["ElapsedSecondsAllTrajectories"]=elapsed_ns/1e9
    m=re.search("@([a-zA-Z0-9\-\.]+):(\w+)\/(\w+)\.git:(\w+)",gitrepo)
    if m:
        website, who, project, version=m.groups()
        info['RepoUrl']='https://{0}/{1}/{2}/commit/{3}'.format(website, who, project, version)
        info['RepoEnglish']="{0}/{1} project at {2}".format(who, project, website)
    s, e, i, r=attrs['Initial Values']
    info['InitialValues']='Susceptible {0} exposed {1} infectious {2}  recovered {3}'.format(
        s, e, i, r)
    root=ET.fromstring(attrs['Options'].tostring().decode("utf-8"))
    commandline_options=dict()
    for option in root.iter("option"):
        name=option.find("name").text
        values=list()
        for v in option.find("values").findall("value"):
            values.append(v.text)
        commandline_options[name]=", ".join(values)

    info['CommandlineOptions']=commandline_options

    first_dset=list(h5f["/trajectory"].keys())[0]
    param_dict=dict()
    model_params=h5f["trajectory/{0}/seirtotal".format(first_dset)].attrs
    for k, v in model_params.items():
        param_dict[k]=v[0]
    info["parameters"]=param_dict

    info['UUID']=attrs['uuid'].tostring().decode("utf-8")

    dsetname=h5i["/images"].attrs["largesttrajectory"].tostring().decode("utf-8")
    init=quickpen.initial_values(h5f, dsetname)
    info['PenCount']=init.shape[0]
    info['TotalAnimals']=s+e+i+r
    info['AnimalsPerPen']=int((s+e+i+r)/init.shape[0])
    return info


def write_report(info, outfile):
    summary_table=["\\begin{tabular}{ll}", "Summary & Value \\\\ \\hline"]
    beta0=float(info["CommandlineOptions"]["beta0"])
    summary_table.append("$R_0$ within pen & {0} \\\\".format(3.585*beta0))
    summary_table.append("Total pens & {0} \\\\".format(info['PenCount']))
    summary_table.append("Total animals & {0} \\\\".format(info['TotalAnimals']))
    summary_table.append("Animals per pen & {0} \\\\".format(info['AnimalsPerPen']))
    pen_rider="yes"
    if "rider" in info["CommandlineOptions"]:
        if int(info["CommandlineOptions"]["rider"])==0:
            pen_rider="no"
    fenceline="yes"
    if "disconnected" in info["CommandlineOptions"]:
        if int(info["CommandlineOptions"]["disconnected"])==1:
            fenceline="no"
    alltoall="no"
    if "otherpens" in info["CommandlineOptions"]:
        if int(info["CommandlineOptions"]["otherpens"])==1:
            alltoall="yes"
    demog="no"
    if "replace" in info["CommandlineOptions"]:
        if int(info["CommandlineOptions"]["replace"])==1:
            demog="yes"
    summary_table.append("Pen rider & {0} \\\\".format(pen_rider))
    summary_table.append("Fenceline contact & {0} \\\\".format(fenceline))
    summary_table.append("Background all-to-all infection & {0} \\\\".format(alltoall))
    summary_table.append("Replacement by pen & {0} \\\\".format(demog))
    summary_table.append("\\end{tabular}")
    summary_table="\n".join(summary_table)
    info["SummaryTable"]=summary_table

    options_table=["\\begin{tabular}{ll}", "Commandline Option & Value \\\\ \\hline"]
    for k, v in sorted(info["CommandlineOptions"].items()):
        options_table.append("{0} & {1} \\\\".format(k, v))
    options_table.append("\\end{tabular}")
    options_table="\n".join(options_table)
    info["CommandlineOptionsTable"]=options_table

    provenance_table=["\\begin{tabular}{ll}", "Trait & Value \\\\ \\hline"]
    provenance_table.append("Compile time & {0} \\\\".format(info["CompileTime"]))
    if "RepoUrl" in info:
        provenance_table.append(
            "Code repository & \\href{{{0}}}{{{1}}} \\\\".format(
                info["RepoUrl"], info["RepoEnglish"]))
    provenance_table.append("Initial values & {0} \\\\".format(info["InitialValues"]))
    provenance_table.append("Unique Tag & {0} \\\\".format(info["UUID"]))
    provenance_table.append("Wall time [s] & {0:7.2f} \\\\".format(info["ElapsedSecondsAllTrajectories"]))
    provenance_table.append("\\end{tabular}")
    info["CodeTraitsTable"]="\n".join(provenance_table)

    params_table=["\\begin{tabular}{ll}", "Internal Parameter & Value \\\\ \\hline"]
    param_names=list(info["parameters"].keys())
    param_names.sort()
    for k in param_names:
        params_table.append("{0} & {1} \\\\".format(k, info["parameters"][k]))
    params_table.append("\\end{tabular}")
    info["ParametersTable"]="\n".join(params_table)

    info["TrajectoryLines"]=include_trajectory_lines()
    info["EndTime"]=include_end_time()
    info["TotalInfected"]=include_total_infected()
    info["TotalInfectedPoints"]=include_total_infected_points()
    info["BinnedPrevalence"]=include_binned()
    info["SmallMultiples"]=include_multiples()
    info["TrajectoryDensityInfected"]=include_trajectory_density()

    text="""\\documentclass{{article}}
\\usepackage{{graphicx}}
\\usepackage[margin=0.5in]{{geometry}}
\\usepackage{{hyperref}}
\\usepackage{{palatino}}
\\begin{{document}}
\\title{{Report on {Title}}}
\\date{{\\today}}
\\author{{Generated with \\texttt{{report.py}}}}
\\maketitle

{SummaryTable}

{CommandlineOptionsTable}

{ParametersTable}

{CodeTraitsTable}

{TrajectoryLines}

{TrajectoryDensityInfected}

{BinnedPrevalence}

{EndTime}

{TotalInfectedPoints}

{TotalInfected}

{SmallMultiples}


\\end{{document}}
""".format(**info)
    f=open(outfile, "w")
    f.write(text)


def single_trajectory_small_multiples(h5f, h5i):
    dsetname=h5i["/images"].attrs["largesttrajectory"].tostring().decode("utf-8")
    #dsetname=quickpen.trajectories(h5f)[0]
    tr,times,obs=quickpen.per_pen_trajectory(h5f, dsetname)
    penplot.small_multiples(tr, times, obs, 4)


def summaries(h5f, h5i):
    #aggregate=quickpen.summary_of_ensemble(h5f, -1)

    total_infected=h5i["/images/totalinfected"][:]
    total_infected.sort()
    end_time=h5i["/images/endtime"][:]
    end_time.sort()

    penplot.total_infected_plot(total_infected)
    penplot.total_infected_count_plot(total_infected)
    penplot.end_time_plot(end_time)
    penplot.survival_susceptible(h5i["/images/infectiontimes"])

    infected_times=h5i["/images/trajectorydensityx"]
    infected_counts=h5i["/images/trajectorydensityy"]
    infected_hist=h5i["/images/trajectorydensity"]
    penplot.trajectory_image(infected_hist, infected_times, infected_counts,
        "Infected")

def include_trajectory_density():
    return include_figure("trajectory_density_Infected.pdf", "0.7",
        "This smooths over all realizations in the ensemble in order "
        +"to create an estimate of the probability distribution "
        +"for finding the system at a given state and time.",
        "fig:trajectorydensityinfected")

def include_multiples():
    return include_figure("multiples.pdf", "0.7",
        "Each subgraph is a separate pen, showing exposed and infected "
        +"over time. This is one sample realization from the file.",
        "fig:smallmultiples");

def include_end_time():
    return include_figure("end_time_hist.pdf", "0.7",
        "Each bar shows a count of how many realizations completed "+
        "at a given time.", "fig:endtimes")

def survival_infected():
    return include_figure("survival_infected_hist.pdf", "0.7",
        "This shows the probability of survival to a time. "+
        "", "fig:endtimes")

def include_total_infected():
    return include_figure("total_infected_hist.pdf", "0.7",
        "Each bar shows the total number of realizations whose "
        "infections fell in the given range", "fig:totalinfected")

def include_total_infected_points():
    return include_figure("total_infected_points.pdf", "0.7",
        "Simple point plot of total infected.", "fig:totalinfectedpoint")

def trajectory_lines(h5f):
    penplot.plot_trajectory_lines(quickpen.FileTrajectories(h5f))

def include_trajectory_lines():
    return include_figure("trajectory_lines_Infected.pdf", "0.7",
        "Each line represents a separate realization from the simulation."+
        "This shows all infected individuals, whether exposed or infectious.",
        "fig:trajlines")

def binned_trajectories(h5f, h5i):
    prevalence=h5i["/images/prevalence"]
    trajectories=quickpen.trajectories(h5f)
    initial=quickpen.initial_values(h5f, trajectories[0])
    initial=np.sum(initial[:,:4], axis=0)
    logger.debug("binned_trajectories initial {0}".format(initial[:]))
    total_individuals=np.sum(initial)
    # only count trajectories that make it to this day
    # denom=np.repeat(np.expand_dims(np.sum(prevalence, axis=1), axis=1),
    #         5, axis=1)
    denom=np.ones(prevalence.shape,
        dtype=np.float)*total_individuals*len(trajectories)
    normed=prevalence/denom
    normed[np.isnan(normed)]=0
    logger.debug("binned_trajectories normed type {0}".format(normed.dtype))
    #logger.debug(normed)
    # binned=np.zeros((200, 5), np.float32)
    # ft=quickpen.FileTrajectories(h5f)
    # cnt=len(ft)
    # for i in range(len(ft)):
    #     total, times=ft[i]
    #     quickpen.add_to_binned_trajectory(binned, total, times)
    # binned/=cnt
    # logger.debug("binned {0}".format(binned))
    penplot.prevalence_by_day(normed)

def include_binned():
    return include_figure("prevalencebyday.pdf", "0.7",
        "Exposed and susceptible counts, averaged over all realizations "
        +"in the ensemble. The small horizontal lines indicate that each "
        +"observation is a daily measurement. Exposed is blue, "
        +"infectious in green.", "fig:prevalencebyday")

def include_figure(fileglob, scale, caption, label):
    snippet=""
    targets=glob.glob(fileglob)
    blanks={"scale" : scale, "caption" : caption}
    for idx, f in enumerate(targets):
        blanks["filename"]=f
        if len(targets)>1:
            blanks["label"]="{0}{1}".format(label, idx)
        else:
            blanks["label"]=label
        snippet=snippet+"""\\begin{{figure}}
\\centerline{{\\includegraphics[scale={scale}]{{{filename}}}}}
\\caption{{{caption}\\label{{{label}}}}}
\\end{{figure}}
""".format(**blanks)

    return snippet


def make_report(filename, h5i):
    f=h5py.File(filename, "r")
    info=dict()
    info.update(metadata(f, h5i))
    fileshort=os.path.basename(filename)
    info["Title"]="{0}--{1}".format(fileshort, info["UUID"][0:5])
    outfile=".".join(os.path.basename(filename).split(".")[:-1]) + ".tex"
    logger.info("writing report to {0}".format(outfile))
    write_report(info, outfile)
    return outfile


def parallel_generate(filename, pyfile, imagefilename, timeout=3.14e7):
    '''
    Timeout in minutes.
    '''
    pyex=sys.executable
    todo=[
      [pyex, pyfile, "--file", filename, "--image", imagefilename, "--multiples"],
      [pyex, pyfile, "--file", filename, "--image", imagefilename, "--summary"],
      [pyex, pyfile, "--file", filename, "--image", imagefilename, "--binned"]
    ]
    processes=list()
    for t in todo:
        processes.append(subprocess.Popen(t))
    interval=2 # seconds
    for check in range(int(timeout*60/interval)):
        done=list()
        for idx, proc in enumerate(processes):
            retcode=proc.poll()
            if retcode!=None:
                done.append(idx)
                logging.info("{0} completed {1}".format(todo[idx][4], retcode))
        done.reverse()
        for d in done:
            processes.pop(d)
            todo.pop(d)
        if len(processes)==0:
            return True
        time.sleep(interval)
    return False


if __name__ == "__main__":
    #logging.basicConfig(level=logging.INFO)
    parser=DefaultArgumentParser(description="Quick look at an H5 file")
    parser.add_argument("--file", dest="file", action="store",
        default="rider.h5", help="data file to read")
    parser.add_argument("--image", dest="image", action="store",
        default="image.h5", help="summary data from file")
    parser.add_function("report", "Build the report")
    parser.add_function("lines", "Write the image of all realizations")
    parser.add_function("binned", "Average daily prevalence")
    parser.add_function("multiples", "Small multiples graph of one realization")
    parser.add_function("summary", "Several summary graphs")
    parser.add_function("generate", "Create graphs and report")

    args=parser.parse_args()

    pyfile=os.path.abspath(__file__)
    filename=os.path.abspath(args.file)
    imagefilename=os.path.abspath(args.image)
    dirname=".".join(os.path.basename(filename).split(".")[:-1])+"d"
    f=h5py.File(filename, "r")
    h5i=h5py.File(imagefilename, "r")
    if args.report:
        texfile=make_report(filename, h5i)
        subprocess.call(["latexmk", "-pdf", texfile])
        sys.exit()
    elif args.generate:
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        os.chdir(dirname)
        parallel_generate(filename, pyfile, imagefilename)
        texfile=make_report(filename, h5i)
        subprocess.call(["latexmk", "-pdf", texfile])
        sys.exit()

    if args.lines:
        trajectory_lines(f)
    elif args.binned:
        binned_trajectories(f, h5i)
    elif args.multiples:
        single_trajectory_small_multiples(f, h5i)
    elif args.summary:
        summaries(f, h5i)
    if not parser.any_function():
        parser.print_help()
