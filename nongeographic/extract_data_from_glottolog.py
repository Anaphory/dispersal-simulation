import xml.etree.ElementTree as ET
import itertools
import pyglottolog

from cldfcatalog import Config
cfg = Config.from_file()
glottolog = pyglottolog.Glottolog(cfg.get_clone("glottolog"))
languoids = {l.id: l for l in glottolog.languoids()}


def all_macroareas(family):
    macroareas = set(family.macroareas)
    for c in family.descendants_from_nodemap(languoids):
        macroareas |= set(c.macroareas)
    return macroareas

families_by_macroarea = {}
for toplevel in glottolog.tree.glob("*"):
    # These are the top-level families, I guess there is a better way to access
    # them, but I don't find it documented or by reading the pyglottolog code.
    glottocode = toplevel.stem
    # A few glottocodes are bookkeeping families, not actual language families
    family = glottolog.languoid(glottocode)
    if family.category == "Pseudo Family":
        continue
    macroareas = all_macroareas(family)
    if len(macroareas) == 1:
        families_by_macroarea.setdefault(
            macroareas.pop().id, []).append(family)
    else:
        # Some language families have a clear belonging to one macroarea, but
        # contain off-branches or errors in other parts of the world. Clear them up.
        macroarea = {
            'japo1237': 'eurasia',
            'indo1319': 'eurasia',
            'aust1307': 'pacific',
            'araw1281': 'southamerica',
            'aust1305': 'eurasia',
            'eski1264': 'northamerica',
            'grea1241': 'eurasia',
            'atla1278': 'africa',
            'maya1287': 'northamerica',
        }.get(family.id)
        if macroarea:
            families_by_macroarea.setdefault(
                macroarea, []).append(family)
        else:
            pass


root = ET.ElementTree(ET.fromstring(
    """<?xml version='1.0' encoding='UTF-8'?>
    <beast beautistatus="" beautitemplate="Standard" namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" version="2.0">
      <map name="taxon">beast.evolution.alignment.Taxon</map>
      <map name="taxonset">beast.evolution.alignment.TaxonSet</map>
      <branchRateModel clock.rate="@clockRate.c:default" id="StrictClockModel.c:default" spec="beast.evolution.branchratemodel.StrictClockModel" />
      <run id="mcmc" spec="MCMC" chainLength="500000000" numInitializationAttempts="1000" sampleFromPrior="false">
        <state id="state" storeEvery="5000">
            <!-- states go here -->
        </state>
        <distribution id="posterior" spec="util.CompoundDistribution">
            <distribution id="prior" spec="util.CompoundDistribution">
            <!-- tree and other priors go here -->
            </distribution>
            <distribution id="likelihood" spec="util.CompoundDistribution">
            </distribution>
        </distribution>
        <logger id="screenlog" logEvery="500">
            <log arg="@posterior" id="ESS.0" spec="util.ESS" />
            <log idref="prior" />
            <log idref="likelihood" />
            <log idref="posterior" />
        </logger>
        <!-- inits go here -->
        <!-- operators go here -->
        <!-- loggers go here -->
      </run>
    </beast>"""))

state = root.find("./run/state")
mcmc = root.find("./run")
prior = root.find("./run/distribution/distribution")

for macroarea, families in families_by_macroarea.items():
    family_sizes = {
        family.id: [
            language.id for language in family.descendants_from_nodemap(
                languoids, level=glottolog.languoid_levels["language"])] +
        ([family.id] if family.level.id == 'language' else [])
        for family in families}

    ET.SubElement(
        state, "parameter",
        id="birthRate.t:{:}".format(macroarea), name="stateNode", value="8e-4")
    ET.SubElement(
        state, "parameter",
        id="deathRate.t:{:}".format(macroarea), name="stateNode", value="5e-1")
    ET.SubElement(
        state, "parameter",
        id="sampling.t:{:}".format(macroarea), name="stateNode", value="2e-1")

    ET.SubElement(
        mcmc, "operator",
        id="deathRateScaler.t:{:}".format(macroarea),
        parameter="@deathRate.t:{:}".format(macroarea),
        scaleFactor="0.5", spec="ScaleOperator", weight="3.0")
    ET.SubElement(
        mcmc, "operator",
        id="birthRateScaler.t:{:}".format(macroarea),
        parameter="@birthRate.t:{:}".format(macroarea),
        scaleFactor="0.5", spec="ScaleOperator", weight="3.0")
    ET.SubElement(
        mcmc, "operator",
        id="samplingScaler.t:{:}".format(macroarea),
        parameter="@sampling.t:{:}".format(macroarea),
        scaleFactor="0.5", spec="ScaleOperator", weight="3.0")

    logger = ET.SubElement(
        mcmc, "logger",
        fileName="{:}.log".format(macroarea),
        logEvery="500")
    ET.SubElement(logger, "log", idref="birthRate.t:{:}".format(macroarea))
    ET.SubElement(logger, "log", idref="deathRate.t:{:}".format(macroarea))
    ET.SubElement(logger, "log", idref="sampling.t:{:}".format(macroarea))

    for family, languages in itertools.chain(
            [(macroarea, ["proto-" + family for family in family_sizes])],
            family_sizes.items()):
        ET.SubElement(
            prior,
            "distribution",
            birthDiffRate="@birthRate.t:{:}".format(macroarea),
            id="BirthDeathModel.t:{:}".format(family),
            relativeDeathRate="@deathRate.t:{:}".format(macroarea),
            sampleProbability="@sampling.t:{:}".format(macroarea),
            spec="beast.evolution.speciation.BirthDeathGernhard08Model",
            tree="@Tree.t:{:}".format(family),
            type="unscaled")

        ET.SubElement(
            ET.SubElement(
                ET.SubElement(
                    mcmc,
                    "init",
                    estimate="false",
                    id="startingTree.t:{:}".format(family),
                    initial="@Tree.t:{:}".format(family),
                    rootHeight="{:}".format(6000.0),
                    spec="beast.evolution.tree.RandomTree",
                    taxonset="@taxa:{:}".format(family)),
                "populationModel",
                spec="ConstantPopulation"),
            "popSize",
            spec="parameter.RealParameter",
            value="1")


        ET.SubElement(
            ET.SubElement(
                ET.SubElement(
                    ET.SubElement(
                        state,
                        "tree",
                        id="Tree.t:{:}".format(family),
                        name="stateNode"),
                    "taxonset",
                    id="taxa:{:}".format(family),
                    spec="TaxonSet"),
                "plate",
                range="{:}".format(','.join(languages)),
                var="language"),
            "taxon",
            id="$(language)")

        logger = ET.SubElement(
            mcmc, "logger",
            fileName="{:}-{:}.nex".format(macroarea, family),
            mode="tree",
            logEvery="500")
        ET.SubElement(logger, "log", idref="Tree.t:{:}".format(family))

        ET.SubElement(
            mcmc, "operator",
            id="SubtreeSlide.t:{:}".format(family),
            markclades="true", spec="SubtreeSlide", tree="@Tree.t:{:}".format(family),
            weight="15.0")
        ET.SubElement(
            mcmc, "operator",
            id="narrow.t:{:}".format(family),
            markclades="true", spec="Exchange", tree="@Tree.t:{:}".format(family),
            weight="15.0")
        ET.SubElement(
            mcmc, "operator",
            id="wide.t:{:}".format(family),
            isNarrow="false", markclades="true", spec="Exchange", tree="@Tree.t:{:}".format(family),
            weight="3.0")
        ET.SubElement(
            mcmc, "operator",
            id="WilsonBalding.t:{:}".format(family),
            markclades="true", spec="WilsonBalding", tree="@Tree.t:{:}".format(family),
            weight="3.0")
        # java.lang.IllegalArgumentException: n must be positive
        # ET.SubElement(
        #     mcmc, "operator",
        #     id="UniformOperator.t:{:}".format(family),
        #     spec="Uniform", tree="@Tree.t:{:}".format(family),
        #     weight="30.0")

        # Deactivate root scaling

        # ET.SubElement(
        #     mcmc, "operator",
        #     id="treeScaler.t:{:}".format(family),
        #     scaleFactor="0.5", spec="ScaleOperator", tree="@Tree.t:{:}".format(family),
        #     weight="3.0")
        # ET.SubElement(
        #     mcmc, "operator",
        #     id="treeRootScaler.t:{:}".format(family),
        #     rootOnly="true", scaleFactor="0.5", spec="ScaleOperator", tree="@Tree.t:{:}".format(family),
        #     weight="3.0")
    break

root.write("beast.xml")
