#include "tools/string_utils.h"

#include "web/web.h"
#include "web/Input.h"
#include "web/emfunctions.h"
// #include "web/d3/d3_init.h"
#include "web/d3/selection.h"
#include "web/d3/svg_shapes.h"

#include "tools/Graph.h"

#include "interaction_networks.h"

namespace UI = emp::web;

// UI::Document doc("emp_base");
Controller c;
D3::PowScale edge_color;
D3::LinearScale edge_width;
D3::LinearScale node_color;
D3::ToolTip edge_tool_tip([](double d) {return D3::FormatFunction(".4f")(d);});
D3::ToolTip node_tool_tip([](std::string d) {return d;});

UI::Document div_controls("div_controls");

void ClearGraph() {
    D3::Select("#lexicase_graph").SelectAll("circle").Remove();
    D3::Select("#lexicase_graph").SelectAll("path").Remove();

    D3::Select("#eco_ea_graph").SelectAll("circle").Remove();
    D3::Select("#eco_ea_graph").SelectAll("path").Remove();

    D3::Select("#sharing_graph").SelectAll("circle").Remove();
    D3::Select("#sharing_graph").SelectAll("path").Remove();

}

void ResetScales() {
    double min_edge = 0;
    double max_edge = 0;

    for (auto w_vec : c.lex_network.GetWeights()) {
        for (double w : w_vec) {
            if (w < min_edge) {
                min_edge = w;
            }
            else if (w > max_edge) {
                max_edge = w;
            }
        }
    }

    for (auto w_vec : c.eco_network.GetWeights()) {
        for (double w : w_vec) {
            if (w < min_edge) {
                min_edge = w;
            }
            else if (w > max_edge) {
                max_edge = w;
            }
        }
    }

    for (auto w_vec : c.share_network.GetWeights()) {
        for (double w : w_vec) {
            if (w < min_edge) {
                min_edge = w;
            }
            else if (w > max_edge) {
                max_edge = w;
            }
        }
    }

    edge_color.Exponent(.25);
    edge_color.SetDomain(emp::array<double, 3>({min_edge, 0, max_edge}));
    edge_color.SetInterpolate("interpolateLab");
    edge_color.SetRange(emp::array<std::string, 3>({"red", "grey", "blue"}));

    edge_width.SetDomain(emp::array<double, 3>({min_edge, 0, max_edge}));
    edge_width.SetRange(emp::array<double, 3>({5, 0, 5}));
}


void DrawGraph(emp::WeightedGraph g, std::string canvas_id, double radius = 150) {
    D3::Selection s = D3::Select(canvas_id);
    D3::Selection defs = D3::Select("#arrow_defs");

    std::function<std::string(std::string)> MakeArrow = [&defs](std::string color) {
        std::string id = "arrow_" + color;
        emp::remove_chars(id, "#,() ");

        defs.Append("svg:marker")
            .SetAttr("id", id)
            .SetAttr("viewBox", "0 -5 10 10")
            .SetAttr("refX", 5) // This sets how far back it sits, kinda
            .SetAttr("refY", 0)
            .SetAttr("markerWidth", 9)
            .SetAttr("markerHeight", 9)
            .SetAttr("orient", "auto")
            .SetAttr("markerUnits", "userSpaceOnUse")
            .Append("svg:path")
            .SetAttr("d", "M0,-5L10,0L0,5")
            .SetStyle("fill", color);
        
        return "url(#" + id + ")";
    };

    double theta = 0;
    double inc = 2*3.14 / g.GetNodes().size();
    for (emp::Graph::Node n : g.GetNodes()) {
        double cx = radius + cos(theta) * radius;
        double cy = radius + sin(theta) * radius;
        theta += inc;
        D3::Selection new_node = s.Append("circle");
        new_node.SetAttr("r", 10).SetAttr("cx", cx).SetAttr("cy", cy).BindToolTipMouseover(node_tool_tip);;
        // std::cout << n.GetLabel() << std::endl;
        EM_ASM_ARGS({js.objects[$0].datum(Pointer_stringify($1));}, new_node.GetID(), n.GetLabel().c_str());
    }

    D3::LineGenerator l;
    auto weights = g.GetWeights();
    for (int i = 0; i < weights.size(); ++i) {
        for (int j = 0; j < weights[i].size(); ++j) {
            // std::cout << weights[i][j] << std::endl;
            if (weights[i][j]) {
                double cxi = radius + cos(i*inc) * radius;
                double cxj = radius + cos(j*inc) * radius;
                double cyi = radius + sin(i*inc) * radius;
                double cyj = radius + sin(j*inc) * radius;
                std::string color = edge_color.ApplyScaleString(weights[i][j]);
                // std::cout << color << std::endl;
                emp::array<emp::array<double, 2>, 2> data({emp::array<double, 2>({cxi, cyi}), emp::array<double, 2>({cxj, cyj})});
                D3::Selection n = s.Append("path");
                n.SetAttr("d", l.Generate(data))
                 .SetStyle("stroke", color)
                 .SetStyle("stroke-width", edge_width.ApplyScale(weights[i][j]))
                 .SetAttr("marker-end", MakeArrow(color))
                 .BindToolTipMouseover(edge_tool_tip);

                EM_ASM_ARGS({js.objects[$0].datum($1);}, n.GetID(), weights[i][j]);
                //  std::cout << MakeArrow(color) << std::endl;
                // s.Append("path").SetAttr("d", "M"+emp::to_string(cxi)+","+emp::to_string(cyi)+"L"+emp::to_string(cxj)+","+emp::to_string(cyj)).SetStyle("stroke", color).SetStyle("stroke-width", weights[i][j]*10);
            }
        }
    }

    s.SetupToolTip(edge_tool_tip);
    s.SetupToolTip(node_tool_tip);

}

// std::string GetEdgeColor(double w) {
//     EM_ASM_ARGS({empjs.objects[$0]}, edge_color.GetID());
// }

int main() {
    
    auto pop_selector = UI::Input([](std::string curr){ 
        c.SetPopSize(emp::from_string<double>(curr));
        // std::cout << c.GetPopSize() << std::endl;
    }, "range", "Population size", "popsize_slider");
    pop_selector.Min(2);
    pop_selector.Max(100);
    pop_selector.Value(10);

    auto trait_selector = UI::Input([](std::string curr){ 
        c.SetNTraits(emp::from_string<double>(curr));
    }, "range", "Number of traits", "ntraits_slider");
    trait_selector.Min(1);
    trait_selector.Max(25);
    trait_selector.Value(1);

    auto sigma_share_selector = UI::Input([](std::string curr){ 
        c.SetSigmaShare(emp::from_string<double>(curr));
        // std::cout << c.GetSigmaShare() << std::endl;
    }, "range", "Sharing threshold", "sigmashare_slider");
    sigma_share_selector.Min(0);
    sigma_share_selector.Max(50);
    sigma_share_selector.Value(8);
    sigma_share_selector.Step(.5);

    auto alpha_selector = UI::Input([](std::string curr){ 
        c.SetAlpha(emp::from_string<double>(curr));
    }, "range", "Alpha", "alpha_slider");
    alpha_selector.Min(0);
    alpha_selector.Max(2);
    alpha_selector.Value(1);
    alpha_selector.Step(.1);

    auto cost_selector = UI::Input([](std::string curr){ 
        c.SetCost(emp::from_string<double>(curr));
    }, "range", "Cost", "cost_slider");
    cost_selector.Min(0);
    cost_selector.Max(10);
    cost_selector.Value(1);
    cost_selector.Step(.1);

    auto cf_selector = UI::Input([](std::string curr){ 
        c.SetCf(emp::from_string<double>(curr));
    }, "range", "Consumption Fraction", "cf_slider");
    cf_selector.Min(0);
    cf_selector.Max(1);
    cf_selector.Value(.0025);
    cf_selector.Step(.0001);

    auto niche_width_selector = UI::Input([](std::string curr){ 
        c.SetNicheWidth(emp::from_string<double>(curr));
    }, "range", "Niche width", "nichewidth_slider");
    niche_width_selector.Min(0);
    niche_width_selector.Max(10);
    niche_width_selector.Value(3);
    niche_width_selector.Step(.5);

    auto max_score_selector = UI::Input([](std::string curr){ 
        c.SetMaxScore(emp::from_string<double>(curr));
    }, "range", "Max consumption", "maxscore_slider");
    max_score_selector.Min(0);
    max_score_selector.Max(10);
    max_score_selector.Value(3);
    max_score_selector.Step(.5);

    div_controls << pop_selector << "<br>";
    div_controls << trait_selector << "<br>";
    div_controls << sigma_share_selector << "<br>";
    div_controls << alpha_selector << "<br>";
    div_controls << cost_selector << "<br>";    
    div_controls << cf_selector << "<br>";    
    div_controls << niche_width_selector << "<br>";    
    div_controls << max_score_selector << "<br>";    
    div_controls << UI::Button( [](){ c.Regenerate(); ClearGraph(); ResetScales(); DrawGraph(c.lex_network, "#lexicase_graph");DrawGraph(c.eco_network, "#eco_ea_graph");DrawGraph(c.share_network, "#sharing_graph"); }, "Redraw", "redraw_button");
    c.Regenerate();
    ResetScales(); DrawGraph(c.lex_network, "#lexicase_graph");DrawGraph(c.eco_network, "#eco_ea_graph");DrawGraph(c.share_network, "#sharing_graph");    
    
    // edge_tool_tip.SetHtml([](int weight){return emp::to_string(weight);});

    emp::web::OnDocumentReady([](){
        EM_ASM(d3.select("#redraw_button").classed("btn btn-primary btn-block", true););
    });

    // doc << div_controls;
    // using pop_t = emp::vector<emp::vector<int>>;
    // pop_t pop = make_pop(r);
    // std::cout << emp::to_string(pop) << std::endl;
    // sharing_fitness(pop, Merge(SigmaShare(3), DEFAULT));
    // calc_competition(pop, sharing_fitness);
}