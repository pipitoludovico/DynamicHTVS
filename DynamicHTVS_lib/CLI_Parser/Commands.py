from DynamicHTVS_lib.CLI_Parser import ArgParser

ap = ArgParser.ArgParser()
commands = ap.argumentParser()

dock = None
restrict = None
md = None
rank = None
parameterize = None
build = None
ligands = None
launchOpenmm = None
calculateScore = None
amber = None
batch = None
consider = None
boxsize = None
poses = None
exclude = None
productionTime = None
purge = None
membraneRestraints = None
selection_ = None

dock = commands.get('dock', None)
restrict = commands.get('restrict', None)

rank = commands.get('rank', False)
parameterize = commands.get('parameterize', False)
build = commands.get('build', False)
ligands = commands.get("ligands", ["smi", "smi"])
launchOpenmm = commands.get("launchOpenmm", False)
calculateScore = commands.get("calculatescore", False)
amber = commands.get("amber", None)
batch = commands.get("batch", None)
consider = commands.get("consider", 20)
boxsize = commands.get("boxsize", 12)
poses = commands.get("poses", 1)
exclude = commands.get("exclude", None)
productionTime = commands.get("prt", 10)
purge = commands.get("purge", False)
membraneRestraints = commands.get("membraneRestraints", None)

if commands.get("selection"):
    selection_ = " ".join(commands.get("selection"))

md = commands.get('md', False)

if md:
    parameterize, build, launchOpenmm, calculateScore = True, True, True, True
