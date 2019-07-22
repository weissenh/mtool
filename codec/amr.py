import re
import sys

from graph import Graph
from smatch.amr import AMR;
from collections import Counter

def amr_lines(fp, alignment):
    id, snt, lines = None, None, [];
    alignment = read_alignment(alignment);
    for line in fp:
        line = line.strip();
        if len(line) == 0:
            if len(lines) > 0:
                i = mapping = None;
                try:
                    i, mapping = next(alignment);
                except:
                    print("amr_lines(): missing alignment for graph #{}."
                          "".format(id), file = sys.stderr);
                    pass;
                yield id, snt, " ".join(lines), \
                    mapping if mapping is not None and i == id else None;
            id, lines = None, []
        else:
            if line.startswith("#"):
                if line.startswith("# ::id"):
                    id = line.split()[2]
                if line.startswith("# ::snt"):
                   snt = line[8:].strip();
            else:
                lines.append(line)
    if len(lines) > 0:
        i = mapping = None;
        try:
            i, mapping = next(alignment);
        except:
            print("amr_lines(): missing alignment for graph #{}."
                  "".format(id), file = sys.stderr);
            pass;
        yield id, snt, " ".join(lines), \
            mapping if mapping is not None and i == id else None;

def read_alignment(stream):
    if stream is None:
        while True: yield None, None;
    else: 
        id = None;
        alignment = dict();
        for line in stream:
            line = line.strip();
            if len(line) == 0:
                yield id, alignment;
                id = None;
                alignment.clear();
            else:
                if line.startswith("#"):
                    if line.startswith("# ::id"):
                        id = line.split()[2];
                else:
                    fields = line.split("\t");
                    if len(fields) == 2:
                        start, end = fields[1].split("-");
                        span = set(range(int(start), int(end) + 1));
                        fields = fields[0].split();
                        if len(fields) > 1 and fields[1].startswith(":"):
                            fields[1] = fields[1][1:];
                            if fields[1] == "wiki": continue;
                        if fields[0] not in alignment: alignment[fields[0]] = bucket = dict();
                        else: bucket = alignment[fields[0]];
                        path = tuple(fields[1:]);
                        if path not in bucket: bucket[path] = can = set();
                        else: can =  bucket[path];
                        can |= span;
        yield id, alignment;

def amr2graph(id, amr, full = False, reify = False, alignment = None):
    graph = Graph(id, flavor = 2, framework = "amr")
    node2id = {}
    i = 0
    for n, v, a in zip(amr.nodes, amr.node_values, amr.attributes):
        j = i
        node2id[n] = j
        top = False;
        for key, val in a:
            if key == "TOP":
                top = True;
        node = graph.add_node(j, label = v, top=top)
        i += 1
        for key, val in a:
            if key != "TOP" \
               and (key not in {"wiki"} or full):
                if val.endswith("¦"):
                    val = val[:-1];
                if reify:
                    graph.add_node(i, label=val)
                    graph.add_edge(j, i, key)
                    i += 1
                else:
                    node.set_property(key, val);

    for src, r in zip(amr.nodes, amr.relations):
        for label, tgt in r:
            normal = None;
            if label == "mod":
                normal = "domain";
            elif label.endswith("-of-of") \
                 or label.endswith("-of") \
                   and label not in {"consist-of" "subset-of"} \
                   and not label.startswith("prep-"):
                normal = label[:-3];
            graph.add_edge(node2id[src], node2id[tgt], label, normal)

    overlay = None;
    if alignment is not None:
        overlay = Graph(id, flavor = 2, framework = "alignment");
        for node in alignment:
            for path, span in alignment[node].items():
                if len(path) == 0:
                    node = overlay.add_node(node2id[node], label = tuple(span));
        for node in alignment:
            i = node2id[node];
            for path, span in alignment[node].items():
                if len(path) == 1:
                    node = overlay.find_node(i);
                    if node is None: node = overlay.add_node(i);
                    node.set_property(path[0], tuple(span));
                elif len(path) > 1:
                    print("amr2graph(): ignoring alignment path {} on node #{} ({})"
                          "".format(path, source, node));

    return graph, overlay;

def convert_amr_id(id):
    m = re.search(r'wsj_([0-9]+)\.([0-9]+)', id)
    if m:
        return "2%04d%03d" % (int(m.group(1)), int(m.group(2)))
    m = re.search(r'lpp_1943\.([0-9]+)', id)
    if m:
        return "1%04d0" % (int(m.group(1)))
    else:
        raise Exception('Could not convert id: %s' % id)

def read(fp, full = False, reify = False,
         text = None, alignment = None, quiet = False):
    n = 0;
    for id, snt, amr_line, mapping in amr_lines(fp, alignment):
        amr = AMR.parse_AMR_line(amr_line)
        if not amr:
            raise Exception("failed to parse #{} ({}); exit."
                            "".format(id, amr_line));
        try:
            if id is not None:
                id = convert_amr_id(id)
            else:
                id = n;
                n += 1;
        except:
            pass
        graph, overlay = amr2graph(id, amr, full, reify, mapping);
        cid = None;
        if text:
            graph.add_input(text, quiet = quiet);
        elif snt:
            graph.add_input(snt, quiet = quiet);
        yield graph, overlay;


def write(g: Graph, stream, mark_props=False):
    # mark_props if True, special marker for properties  _prop  or something
    # todo revertible? see https://github.com/cfmrp/mtool/issues/35
    # e.g. :wiki is lost during read in, node names kinda arbitrary for amr?
    # Note: no alignments in output, requires unique top node
    # todo are there any nodes/edges that need special treatment?
    # todo: print other field info of nodes and edges?
    if g.flavor != 2:  # alignments cannot be considered!
        raise ValueError("AMR output cannot represent alignments! Graph ID = "
                         + str(g.id))
    g.prepare_4_dfsearch()  # invert edges where necessary
    # this preparation can raise exception (not 1 top node, disconnected)
    # todo: really change graph or have a copy of the graph?
    tops = [node for node in g.nodes if node.is_top]
    if len(tops) != 1:  # todo: should we allow empty graph?
        raise ValueError("Printing graph to penman style requires unique top "
                         "node!")
    topnode = tops[0]
    visited_edges = set()
    visited_node2label = dict()
    firstcharcntr = Counter()
    nodestr, visited_node2label, firstcharcntr = get_node_string(node=topnode,
         visted_node2label=visited_node2label, charcntr=firstcharcntr,
         mark_props=mark_props)
    outstr = "("
    outstr += nodestr  # w / work-01
    intend = " " * 5
    outstr += get_string_for_outgoing_edges(graph=g, node=topnode,
        intend=intend, visited_node2label=visited_node2label,
        charcntr=firstcharcntr, visited_edges=visited_edges,
        mark_props=mark_props)
    print(outstr + ")\n", file=stream)
    # todo: this shouldn't happen (since prepare4dfsearch), change to assertion
    if not (visited_node2label.keys() == set(g.nodes)):
        # unvisited_nodes = set(g.nodes) - visited_node2label.keys()
        raise Warning("Some graph nodes weren't covered! Graph ID = " +
                      str(g.id))
    if not (visited_edges == g.edges):
        # unvisited_edges = g.edges - visited_edges
        raise Warning("Some graph edges weren't covered! Graph ID = " +
                      str(g.id))
    return


def get_node_string(node, visted_node2label: dict, charcntr: Counter,
                    mark_props: bool) -> tuple:
    # (1)   w  (seen before)  or
    # (2)   w / want-01  (not seen)  or
    # (3)   n / name :op1 "Pierre" :op2 "Vinken"
    # todo: make sure node.label doesn't contain  :/()  ?
    # todo: make charcntr stored permanently in function?
    nodestr = ""
    if node in visted_node2label:  # seen before, get label from dictionary
        shortlabel = visted_node2label[node]
        nodestr = shortlabel  # e.g.  w
    else:  # not seen before
        # 1. get   w / want-01
        nodelabel = node.label
        assert (len(nodelabel) > 0)
        firstletter = nodelabel[0].lower()  # todo: check whether allowed char?
        shortlabel = firstletter  # e.g.   w
        if charcntr[firstletter] != 0:
            shortlabel += str(charcntr[firstletter])  # e.g.  w1
        visted_node2label[node] = shortlabel
        charcntr[firstletter] += 1
        nodestr += shortlabel + " / " + nodelabel    # e.g.  w / want-01
        # 2. further add property, value pairs like   :op1 "Pierre"
        assert ((node.properties is None and node.values is None) or
                (len(node.properties) == len(node.values)))
        if node.properties is not None and node.values is not None:
            for prop, value in zip(node.properties, node.values):
                # todo when to add " " around value?
                # :op1 "Pierre"  but  :polarity -  and  :day 29
                # todo: make sure prop val don't contain reserved char? \"
                if mark_props:
                    prop = prop + "-prop"
                nodestr += " :" + prop + " "
                if value.isdigit() or value == "-":  # :polarity -   :day 29
                    nodestr += value
                else:  # :op1: "Pierre"
                    nodestr += "\"" + value + "\""
    return nodestr, visted_node2label, charcntr


def get_string_for_outgoing_edges(graph: Graph, node, intend: str,
                                  visited_node2label: dict, charcntr: Counter,
                                  visited_edges: set, mark_props: bool) -> str:
    # recursive function
    assert (node in visited_node2label)
    # todo: one line or multiline? if multiline which intend?
    s = ""
    # if node.isleaf():  # base case: no outgoing edges - nothing to print
    #     # not necessary, since loop zero times executed, simply returns
    #     return s
    # sorted() just to have canonical order
    for edge in sorted(node.outgoing_edges):  # recursive case
        visited_edges.add(edge)
        assert (edge.src == node.id)
        target = graph.nodes[edge.tgt]  # todo check if IndexError!
        s += "\n" + intend + " :" + edge.lab  # e.g.  :arg0
        # todo: replace /:() in label name
        unvisited_target = target not in visited_node2label
        s += " "
        nodestr, visited_node2label, firstcharcntr = get_node_string(
            node=target, visted_node2label=visited_node2label,
            charcntr=charcntr, mark_props=mark_props)
        if unvisited_target:
            s += "("
        s += nodestr   # e.g.  w  or  w / want-01  or  n /name :op1 "Hans"
        if unvisited_target:
            # if has outgoing edges, proceed with those, then close bracket
            s += get_string_for_outgoing_edges(
                graph=graph, node=target, intend=intend + " " * 6,
                visited_node2label=visited_node2label, charcntr=charcntr,
                visited_edges=visited_edges, mark_props=mark_props)
            s += ")"  # + "\n"
    return s
