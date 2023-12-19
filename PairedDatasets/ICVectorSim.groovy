@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')
@Grab(group='org.codehaus.gpars', module='gpars', version='1.1.0')

import java.net.*
import org.openrdf.model.vocabulary.*
import slib.sglib.io.loader.*
import slib.sml.sm.core.metrics.ic.utils.*
import slib.sml.sm.core.utils.*
import slib.sglib.io.loader.bio.obo.*
import org.openrdf.model.URI
import slib.graph.algo.extraction.rvf.instances.*
import slib.sglib.algo.graph.utils.*
import slib.utils.impl.Timer
import slib.graph.algo.extraction.utils.*
import slib.graph.model.graph.*
import slib.graph.model.repo.*
import slib.graph.model.impl.graph.memory.*
import slib.sml.sm.core.engine.*
import slib.graph.io.conf.*
import slib.graph.model.impl.graph.elements.*
import slib.graph.algo.extraction.rvf.instances.impl.*
import slib.graph.model.impl.repo.*
import slib.graph.io.util.*
import slib.graph.io.loader.*
import groovyx.gpars.GParsPool

System.setProperty("jdk.xml.entityExpsansionLimit", "0");
System.setProperty("jdk.xml.totalEntitySizeLimit", "0");

def factory = URIFactoryMemory.getSingleton()
def annotationsPath = args[0] // "data/gene_annotations.tab";
def resSimPath = args[1] // "data/sim.txt";

println(annotationsPath + "\t" + resSimPath)

class Gene {

  int id
  Set annotations

  public Gene(id, annotations) {
    setId(id)
    setAnnotations(annotations)
  }

  void addAnnotation(annotation) {
    annotations.add(annotation);
  }

  def getAnnotations() {
    annotations
  }

}


def getOntology = {

  URI graph_uri = factory.getURI("http://purl.obolibrary.org/obo/")
  G graph = new GraphMemory(graph_uri)

  // Load OBO file to graph "go.obo"
  GDataConf goConf = new GDataConf(GFormat.RDF_XML, "data/go-basic.owl")
  GraphLoaderGeneric.populate(goConf, graph)

  // Add virtual root for 3 subontologies__________________________________
  URI virtualRoot = factory.getURI("http://purl.obolibrary.org/obo/virtualRoot")
  graph.addV(virtualRoot)

  new File(annotationsPath).splitEachLine('\t') { items ->
    String geneId = items[0];
    for (int i = 1; i < items.size(); i++) {
      URI idURI = factory.getURI("http://purl.obolibrary.org/obo/" + geneId);
      String go_id = items[i];
      URI goURI = factory.getURI("http://purl.obolibrary.org/obo/" + go_id);
      Edge e = new Edge(idURI, RDF.TYPE, goURI);
      graph.addE(e);
    }
  }

  GAction rooting = new GAction(GActionType.REROOTING)
  rooting.addParameter("root_uri", virtualRoot.stringValue())
  GraphActionExecutor.applyAction(factory, rooting, graph)
  return graph
}

def getURIfromName = { name ->
    def id = name.split('\\:')
    name = id[0] + "_" + id[1]
    return factory.getURI("http://purl.obolibrary.org/obo/$name")
}

def getGenes = {
    def genes = []
    new File(annotationsPath).splitEachLine('\t') { items ->
        def annots = new Gene(0, new LinkedHashSet())
        for (int j = 1; j < items.size(); j++) {
            annots.addAnnotation(getURIfromName(items[j]))
        }
        genes.push(annots)
    }
    return genes
}

graph = getOntology()
genes = getGenes()

def sim_id = 0 //this.args[0].toInteger()

SM_Engine engine = new SM_Engine(graph)

String[] flags = [
  SMConstants.FLAG_SIM_GROUPWISE_BMA,
  SMConstants.FLAG_SIM_GROUPWISE_BMM,
  SMConstants.FLAG_SIM_GROUPWISE_MAX,
  SMConstants.FLAG_SIM_GROUPWISE_MIN,
  SMConstants.FLAG_SIM_GROUPWISE_MAX_NORMALIZED_GOSIM
]

String[] pairFlags = [
  SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995,
  SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_SCHLICKER_2006,
  SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998,
  SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_JIANG_CONRATH_1997_NORM
]

ICconf icConf = new IC_Conf_Corpus("ResnikIC", SMConstants.FLAG_IC_ANNOT_RESNIK_1995_NORMALIZED);
String flagGroupwise = flags[sim_id.intdiv(pairFlags.size())];
String flagPairwise = pairFlags[sim_id % pairFlags.size()];
SMconf smConfGroupwise = new SMconf(flagGroupwise);
SMconf smConfPairwise = new SMconf(flagPairwise);
smConfPairwise.setICconf(icConf);

Map<URI, Double> ics = engine.getIC_results(icConf);

double globalAvgIc = 0
int totalAnnotations = 0

def outputFile = new File(resSimPath)
def out = new BufferedWriter(new FileWriter(outputFile))

for(gene: genes) {
    List<String> icValuesForGene = [] // List to store IC values for the gene

    for (go_uri: gene.getAnnotations()) {
        double icValue = ics.getOrDefault(go_uri, 0.0)
        icValuesForGene.add(icValue.toString())
    }

    // Join the IC values with a tab and write them to the file
    String icValuesLine = icValuesForGene.join('\t') + '\n'
    out.write(icValuesLine)
}

out.flush()
out.close()
