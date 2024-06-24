@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')
@Grab(group='org.codehaus.gpars', module='gpars', version='1.1.0')

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
import slib.graph.algo.accessor.GraphAccessor;

System.setProperty("jdk.xml.entityExpsansionLimit", "0");
System.setProperty("jdk.xml.totalEntitySizeLimit", "0");

def factory = URIFactoryMemory.getSingleton()
def annotationsPath = args[0] // "data/gene_annotations.gaf";
def resICPath = args[1] // "data/ic.txt";

println(annotationsPath + "\t" + resICPath)

def getURIfromName = { name ->
    def id = name.split('\\:')
    name = id[0] + "_" + id[1]
    return factory.getURI("http://purl.obolibrary.org/obo/$name")
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
    String geneId = items[1];
    URI idURI = factory.getURI("http://purl.obolibrary.org/obo/" + geneId);
    String go_id = items[4];
    URI goURI = getURIfromName(go_id);
    Edge e = new Edge(idURI, RDF.TYPE, goURI);
    graph.addE(e);
  }

  GAction rooting = new GAction(GActionType.REROOTING)
  rooting.addParameter("root_uri", virtualRoot.stringValue())
  GraphActionExecutor.applyAction(factory, rooting, graph)
  return graph
}


graph = getOntology()

SM_Engine engine = new SM_Engine(graph)

ICconf icConf = new IC_Conf_Corpus("ResnikIC", SMConstants.FLAG_IC_ANNOT_RESNIK_1995_NORMALIZED);

Set<URI> goTerms = GraphAccessor.getClasses(graph)
System.out.println("GO terms : " + goTerms.size())


def out = new PrintWriter(new BufferedWriter(
  new FileWriter(resICPath)))
for (URI goTerm: goTerms) {
    def ic = engine.getIC(icConf, goTerm)
    out.println(goTerm.toString() + "\t" + ic)
}
out.flush()
out.close()
