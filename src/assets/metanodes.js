function draw_bounding_box(title,elems_pos){
    var l = elems_pos[0]['x']
    var r = elems_pos[0]['x']
    var t = elems_pos[0]['y']
    var b = elems_pos[0]['y']
    for (var i =1;i<elems_pos.length;++i){
        if(elems_pos[i]["x"]<l)
            l= elems_pos[i]["x"];
        if(elems_pos[i]["x"]>r)
            r= elems_pos[i]["x"];
        if(elems_pos[i]["y"]<b)
            b= elems_pos[i]["y"];
        if(elems_pos[i]["y"]>t)
            t= elems_pos[i]["y"];
    }
    canvas = document.querySelector("#detail_graph canvas:last-child");
    var ctx = canvas.getContext("2d");
    ctx.strokeRect(l,b,r-l,t-b);
    console.log(elems_pos,l,r,b,t,canvas);
    
}
function export_detail_data(elemid,cy){
    if (cy==undefined){
        cy = document.getElementById(elemid)['_cyreg']["cy"];
    }
    var node_jsons = cy.nodes().jsons();
    var edge_jsons = cy.edges().jsons();
    var output_json = {"gene_groups":[],"pathways":[]}
    var genes = {}
    var zones = {}
    var pathways = {}
    node_jsons.forEach(node => {
        if(node["data"]["Signatures"]!=undefined){
            let zone = node["data"]["Signatures"].join(";");
            if(zones[zone]===undefined){
                zones[zone]={"id":zone,"genes":[node["data"]["id"]]}
            }else{
                zones[zone].genes.push(node["data"]["id"]);
            }
            genes[node["data"]["id"]]=node["data"]["Signatures"];
        }else{
            if(node["data"]["is_pathway"]){
                pathways[node["data"]["id"]]={"id":node["data"]["id"],"name":node["data"]["name"],"ReactomeLink":node["data"]["ReactomeLink"],connected_genes:[]}
            }
        }
    });
    edge_jsons.forEach(edge=> {
        let src = genes[edge["data"]["source"]];
        let tgt = pathways[edge["data"]["target"]];
        console.log(src,tgt,edge["data"]["source"],edge["data"]["target"])
        tgt.connected_genes.push({"id":edge["data"]["source"],"signatures":src});
    });
    for (const p in zones) {
        const element = zones[p];
        output_json["gene_groups"].push(element);            
    }
    for (const p in pathways) {
        const element = pathways[p];
        output_json["pathways"].push(element);            
    }
    let a = document.createElement("a");
    a.download= elemid == "detail_graph"?"multi_signature_data.json":"mono_signature_data.json";
    a.href="data:text/json;charset=utf-8,"+encodeURIComponent(JSON.stringify(output_json));
    a.click();

}
function export_overview_data(){
    const cy = document.getElementById("overview_graph")['_cyreg']["cy"];
    var node_jsons = cy.nodes().jsons();
    var edge_jsons = cy.edges().jsons();
    var output_json = {"Signatures":[],"Intersections":[]}
    node_jsons.forEach(node => {
        if(node["data"]["Cancer"] !==undefined){
            output_json["Signatures"].push({
                "Cancer":node["data"]["Cancer"],
                "Comparison":node["data"]["Comparison"],
                "Filter":node["data"]["Filter"],
                "id":node["data"]["id"],
                "gProfiler":node["data"]["gProfiler"],
                "genes":node["data"]["Signature"]
            });
        }
    });
    edge_jsons.forEach(edge=> {
        if(edge["data"]["type"]=="signature"){
            let commonGenes = {};
            for (let i = 0; i < edge["data"]["elems"].length; i++) {
                const element = edge["data"]["elems"][i];
                const symbol = edge["data"]["symbols"][i].props.children;
                commonGenes[element]=symbol;
                
            }
            output_json["Intersections"].push({
                "Signatures":[edge["data"]["source"],edge["data"]["target"]],
                "commonGenes":commonGenes
            })
        }
    });

    let a = document.createElement("a");
    a.download= "overview_data.json";
    a.href="data:text/json;charset=utf-8,"+encodeURIComponent(JSON.stringify(output_json,space="  "));
    a.click();

}
function display_one_legend(legend_data,canvas,height,width){
    canvas.height=height;
    canvas.width=width;
    let ctx =canvas.getContext("2d");
    ctx.strokeStyle = "black";
    var items =0;
    var spaces=1;
    const itemSize=15;
    const space=5;

    ctx.textBaseline = "middle";

    for(let legendItem in legend_data.genes ){
        ctx.fillStyle = legend_data.genes[legendItem];
        ctx.beginPath();

        ctx.arc(20,itemSize*(items+0.5)+spaces*space,itemSize/2,0,2*Math.PI);
        ctx.stroke();
        ctx.fill();
        ctx.fillStyle = "black";
        
        ctx.fillText(legendItem,40,itemSize*(items+0.5)+spaces*space);
        items+=1;
        spaces+=1;
    }
    for(let legendItem in legend_data.pathway){
        ctx.fillStyle = legend_data.pathway[legendItem].color;
        switch(legend_data.pathway[legendItem].shape){
            case "triangle":
                ctx.beginPath();
                ctx.moveTo(20,spaces*space+(items)*itemSize);
                ctx.lineTo(20+itemSize*Math.sqrt(3)/4,spaces*space+items*itemSize+itemSize*Math.sqrt(3)/2);
                ctx.lineTo(20-itemSize*Math.sqrt(3)/4,spaces*space+items*itemSize+itemSize*Math.sqrt(3)/2);
                ctx.stroke();
                ctx.fill();
                break;
            default:
                ctx.fillRect(10,20*(i)-10,20,20);
        }
        ctx.fillStyle = "black";
        ctx.fillText(legendItem,40,spaces*space+(items+0.5)*itemSize);
        items+=1;
        
        spaces+=1;
    }
    for(let legendItem in legend_data.diseases){
        ctx.fillStyle = legend_data.diseases[legendItem];

        ctx.fillRect(20-itemSize/2,spaces*space+(items)*itemSize,itemSize,itemSize);
        ctx.fillStyle = "black";
        ctx.fillText(legendItem,40,spaces*space+(items+0.5)*itemSize);

        items+=1;
        spaces+=1;

    }
    for(let legendItem in legend_data.comparisons){

        ctx.fillText(legendItem,40,spaces*space+(items+0.5)*itemSize);
        let img = document.createElement("img");
        img.src=legend_data.comparisons[legendItem]["background-image"];
        let n=2;
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                ctx.drawImage(img,20-itemSize/2+j*itemSize/n,spaces*space+(items)*itemSize+i*itemSize/n,itemSize/n,itemSize/n );

            }
        }

        items+=1;
        spaces+=1;

    }


}
function display_legend(multilegend_data,mono_legend_data){
    const multi_canvas = document.querySelector(".legend_canvas #multi_canvas");
    const mono_canvas = document.querySelector(".legend_canvas #mono_canvas");
    const height = Math.max(multi_canvas.clientHeight,mono_canvas.clientHeight);
    const width = Math.max(multi_canvas.clientWidth,mono_canvas.clientWidth);
    // console.log(legend_data);
    // console.log(dash_clientside.callback_context.triggered_id);
    // let canvas = document.querySelector(".legend_canvas #"+dash_clientside.callback_context.triggered_id.split("_")[0]+"_canvas");
    display_one_legend(multilegend_data,multi_canvas,height,width);
    display_one_legend(mono_legend_data,mono_canvas,height,width);
    return dash_clientside.no_update;
}
function tooltip(mouseoverNodeData,mouseoverEdgeData,fake_graph_size,elements,extent,stylesheets,pos_store_data,cur_children,cur_show,cur_bbox,cur_direction){
    if(mouseoverNodeData==null && mouseoverEdgeData==null){
    
    //TODO keep cur if close else reset
        //return [cur_children,cur_show,cur_bbox,cur_direction];
        return [[],false,{},"right",""];
    }else{
        try{
            elem_id = mouseoverEdgeData != null?  mouseoverEdgeData["id"]: mouseoverNodeData["id"];
            const cy = document.getElementById(dash_clientside.callback_context.triggered_id)['_cyreg']["cy"];
            const cy_elem = cy.elements("#" + elem_id).toArray();
            cy.nodes("#" + elem_id).closedNeighborhood().flashClass("testFlash",1000);
            cy.edges("#" + elem_id).sources().flashClass("testFlash",1000);
            cy.edges("#" + elem_id).targets().flashClass("testFlash",1000);

            var direction="left";
            var elem = undefined;
            console.log(cy_elem)
            var content = "";
            let triggered_id =dash_clientside.callback_context.triggered_id;
            if(cy_elem[0].isEdge()){
                if(dash_clientside.callback_context.triggered_id=="overview_graph"){
                    content=[{'type': 'H6', 'namespace': 'dash_html_components', 'props': {'children': "Intersection of "+cy_elem[0].data("source")+" and "+cy_elem[0].data("target")}}].concat(cy_elem[0].data("symbols"))
                }else{
                    content = "edge_tooltip";
                    signatures = cy.elements("#"+cy_elem[0].data('source')).data("Signatures");

                    content = [];
                    content.push( {'type': 'H6', 'namespace': 'dash_html_components', 'props': {'children': "Gene Signatures"}})
                    signature_li = []
                    for (var s of signatures){
                        signature_li.push( {'type': 'Li', 'namespace': 'dash_html_components', 'props': {'children': s}})
                    }
                    content.push(
                        {'type': 'Ul', 'namespace': 'dash_html_components', 'props': {'children': signature_li}}
                    )
                    content.push( {'type': 'H6', 'namespace': 'dash_html_components', 'props': {'children': "Gene symbols : "+cy.elements("#"+cy_elem[0].data('source')).data("tooltip_content")[0].props["children"]}})

                    content.push( {'type': 'H6', 'namespace': 'dash_html_components', 'props': {'children': "Pathway : "+cy.elements("#"+cy_elem[0].data('target')).data("tooltip_content")[0].props["children"]}})
                }
                rmp = cy_elem[0].renderedMidpoint()
                mp = cy_elem[0].midpoint()
                setTimeout(() => document.getElementById(triggered_id+"_tooltip").classList.add("pouet"),500);
                return [content,true,{"x0":rmp['x'],"y0":rmp['y'],"x1":rmp['x'],"y1":rmp['y']},direction,""];
            }else{
                // if( mouseoverNodeData["id"] in pos_store_data){
                //     elem=pos_store_data[mouseoverNodeData["id"]];
                // }else{
                //     for(var i of elements){
                //         if( i["data"]["id"] == mouseoverNodeData["id"]){
                //             elem = i;
                //             break;
                //         }
                //     }
                // }
                if(dash_clientside.callback_context.triggered_id!="overview_graph"){
                    var detail_graph = document.getElementById(dash_clientside.callback_context.triggered_id);
                    var detail_resize_span = document.getElementById(dash_clientside.callback_context.triggered_id+'_span');
                    var x = 0;
                    var y=0;
                    var width=0;
                    var height=0;
                    var wratio=detail_graph.clientWidth/(extent["x2"]-extent["x1"]);
                    var hratio=detail_graph.clientHeight/(extent["y2"]-extent["y1"]);
                    x = cy_elem[0].renderedPosition("x")//+document.getElementById("detail_resize_span").clientWidth
                    y = cy_elem[0].renderedPosition("y")
                    direction = (x/detail_graph.clientWidth>0.5)?"left":"right";
                    if( "weight" in mouseoverNodeData){
                        width = mouseoverNodeData["weight"];
                        height = mouseoverNodeData["weight"];
                    }else{
                        var polygon = [];
                        var zeros = [];
                        var y_crosses = [];
                        for(var i of stylesheets){
                            if(i["selector"] == "node#"+mouseoverNodeData["id"]){
                                width = parseFloat(i["style"]["width"]);
                                height = parseFloat(i["style"]["height"]);
                                var pos = i["style"]["shapePolygonPoints"].split(" ").map(parseFloat);
                                for(var p =0;p<pos.length;p+=2){
                                    polygon.push([pos[p],pos[p+1]]);
                                    if(pos[p+1]==0.0){
                                        zeros.push(p/2);
                                        y_crosses.push(pos[p+1]);
                                    }
                                }

                                break;
                            }
                        }
                        x += (polygon[0][0]/2)*width*wratio;
                        y += (polygon[0][1]/2)*height*hratio;
                        direction=polygon[0][0]<0.0?"left":"right";
                        if(Math.abs(polygon[0][1])*height>Math.abs(polygon[0][0])*width){
                            direction=polygon[0][1]<0.0?"top":"bottom";
                        }
                        width=0
                        height=0
                    }
                    setTimeout(() => document.getElementById(triggered_id+"_tooltip").classList.add("pouet"),500);
                    return [mouseoverNodeData["tooltip_content"],true,{"x0":x-width*wratio/2,"y0":y-height*hratio/2,"x1":x+width*wratio/2,"y1":y+height*hratio/2},direction,""];
                }else{
                    console.log(cy_elem[0].renderedBoundingBox());
                    console.log(cy_elem[0]);
                    let bb = cy_elem[0].renderedBoundingBox({'includeLabels':false});
                    setTimeout(() => document.getElementById(triggered_id+"_tooltip").classList.add("pouet"),500);
                    return [mouseoverNodeData["tooltip_content"],true,{"x0":bb.x1,"y0":bb.y1,"x1":bb.x2,"y1":bb.y2},"right",""];
                }  
            }
        }catch(error){
            console.log(error);
            return [[],false,{},"right",""];

        }
    }
}
function layout_overview(elems){
    function idealEdgeLength(edge){
            if(edge.data("type")=="fake" || edge.data("type")=="signature"){
                return 4;
            }else{
                return 128 ;
            }
    }
    function idealEdgeLengthd3(edge){
        if(edge.data("type")=="fake" || edge.data("type")=="signature"){
            return 2;
        }else{
            return 1;
        }
}
    // document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true,"animate":false,"idealEdgeLength":idealEdgeLength}).run();
    // document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"d3-force","animate":false,"linkStrength":idealEdgeLengthd3}).run();
    document.getElementById("overview_graph")['_cyreg']["cy"].layout({"name":"cola","nodeDimensionsIncludeLabels":true,"animate":false,"edgeLength":idealEdgeLength}).run();

    return dash_clientside.no_update;

}
function highlight_pathway_neighbourhood(nodes){
    const cy = document.getElementById("detail_graph")['_cyreg']["cy"];
    for(var p of nodes){
        if (p["is_pathways"]){
            closedNeighborhood = cy.nodes("#"+p["id"]).closedNeighborhood();

        }
    }
}

function use_elem_weight(ele){
    console.log(ele,ele.data("weight"));
    return ele.data("weight");
}

function box_plots_stats(relayoutData,figure,stats_data,do_stats_data){
        setTimeout(function(){
            console.log(stats_data);
            shapes = [];
            const plot = document.querySelector("#activation_boxplot div.js-plotly-plot");
            const plot_max = get_y_max(plot);
            console.log("plot_max",plot_max);
            var  y_max =0;
            let lines = []
            let xs = []
            for(var i =0;i<stats_data["stats"].length;i++){
                let {x0,x1} = stat_shape_x(plot,stats_data["stats"][i][0],stats_data["stats"][i][1],stats_data["stats"][i][2],stats_data["stats"][i][3]);

                if(x0!=undefined && x1!=undefined){
                    xs.push({x0:x0,x1:x1,label:stats_data["stats"][i][4]})
                }
            }
            xs.sort((a,b) => (a.x1-a.x0)>(b.x1-b.x0)?1:(b.x0-a.x0))
            for(var i =0;i<xs.length;i++){
                let {x0,x1,label} = xs[i];
                let lineIndex=undefined;
                for(var l =0;l<lines.length&&(lineIndex==undefined);l++){
                    let found = false;
                    for(var l2=0;l2<lines[l].length&&!found;l2++){
                        if(x0<=lines[l][l2][1]&&x1>=lines[l][l2][0]){
                            found=true;
                        }
                    }
                    if(!found)
                        lineIndex=l;
                }
                if(lineIndex==undefined){
                    lineIndex=lines.length;
                    lines.push([[x0,x1]]);
                }else{
                    lines[lineIndex].push([x0,x1]);

                }
                let y = create_shape(plot_max,0.02,lineIndex,x0,x1,label,shapes);
                // let y = createStatLink(stats_data["stats"][i][0],stats_data["stats"][i][1],stats_data["stats"][i][2],stats_data["stats"][i][3],0.02,i,stats_data["stats"][i][4],shapes);
                if (y_max<y){
                        y_max=y;
                }
            }
            
            if(shapes.length>0){
                // var y_max = shapes[0].y1;
                // for(var i =1;i<shapes.length;i++){
                //     if(y_max<shapes[i].y1)
                //         y_max = shapes[i].y1;
                // }
                var new_layout = {"shapes":shapes};
                if(plot.layout.yaxis.range[1]<y_max){
                    new_layout["yaxis"]={"range":[plot.layout.yaxis.range[0],y_max]};
                }
                Plotly.relayout(plot,new_layout);
            }else{
                Plotly.relayout(plot,{"shapes":[]});

            }
        },1);
}

function tapMultiSignPathway(pathway_id){
    const cy = document.getElementById("detail_graph")['_cyreg']["cy"];
    cy.edges(".tapped").removeClass("tapped")
    if(pathway_id!=undefined){
        cy.nodes("#"+pathway_id).connectedEdges().addClass("tapped")
    }

}
async function draw_offscreen(graph_id,filename,legend_canvas){
    const cy = document.getElementById(graph_id)['_cyreg']["cy"];

    let scale=Math.max(Math.ceil(3000/cy.width()),Math.ceil(3000/cy.height()));
    let blob_promise = cy.png({full:true,scale:scale,output:"blob-promise"});
    let oc = new OffscreenCanvas(cy.width()*scale,cy.height()*scale);
    
    let outCtx = oc.getContext("2d");
    
    let {x1,x2,y1,y2,w,h} = cy.elements().renderedBoundingBox()

    
    let promise = blob_promise.then((blob)=> createImageBitmap(blob));
    if( legend_canvas!==undefined){
        promise = promise.then((main_blob) => {
            let inLegendCtx = legend_canvas.getContext("2d");
            return createImageBitmap(inLegendCtx.getImageData(0,0,legend_canvas.width,legend_canvas.height)).then((legend_bitmap) => {
                // let scale = Math.floor(oc.height/(4*legend_bitmap.height))
                outCtx.drawImage(legend_bitmap,0,0,legend_bitmap.width,legend_bitmap.height,0,0,legend_bitmap.width*scale,legend_bitmap.height*scale);
                return main_blob;
            });
        })
    }
    await promise.then(imageBitMap => {
        outCtx.drawImage(imageBitMap,0,0,imageBitMap.width,imageBitMap.height,oc.width/2-w*scale/2,oc.height/2-h*scale/2,w*scale,h*scale);

        oc.convertToBlob({"type":"image/png"}).then((blob) => {
            let img =URL.createObjectURL(blob);

            let a = document.createElement("a");
            a.download= filename;
            a.href=img;
            a.click();

        });
    })

}
function on_box_plot_click(clickData,relayoutData,figure,stats_data){
    const plot = document.querySelector("#activation_boxplot div.js-plotly-plot");

    if(clickData!=undefined && Object.hasOwn(clickData,"points") && clickData.points.length>0){
        document.querySelectorAll("#activation_boxplot svg.main-svg g.cartesianlayer g.subplot").forEach(n => n.classList.add("masked"));

        var traceIndex = clickData.points[0].curveNumber;
        var x = clickData.points[0].x;
        data_indices = stats_data["stats"]["data_indices"][stats_data["stats"]["curve_numbers"][traceIndex]+"_"+x];
        const y_axis = stats_data["stats"]["shapes"][data_indices[0]].yref.split(" ")[0];
        let subplot_class = y_axis.replace("y","x")+ y_axis;
        console.log(subplot_class,"#activation_boxplot svg.main-svg g.cartesianlayer g.subplot."+subplot_class);
        document.querySelector("#activation_boxplot svg.main-svg g.cartesianlayer g.subplot."+subplot_class).classList.remove("masked");
        Plotly.relayout(plot,Object.assign({},plot.layout,{"shapes":data_indices.map(i => stats_data["stats"]["shapes"][i])}));
    //     let legendGroup = document.querySelector("#activation_boxplot svg.main-svg g.legend g.scrollbox g.groups:nth-of-type("+(traceIndex+1)+")");
    //     let traceName = legendGroup.__data__[0][0].trace.name;
    //     let filtered_stats_data = [];
    //     stats_data["stats"].forEach(element => {
    //         if((element[0]==traceName && element[1]==x)||(element[2]==traceName && element[3]==x)){
    //             filtered_stats_data.push(element);
    //         }
    //     });
    //     box_plots_stats(relayoutData,figure,{"stats":filtered_stats_data});
    //     console.log(traceIndex,x,traceName,filtered_stats_data);
    }else{
        document.querySelectorAll("#activation_boxplot svg.main-svg g.cartesianlayer g.subplot.masked").forEach(n => n.classList.remove("masked"));
        Plotly.relayout(plot,Object.assign({},plot.layout,{"shapes":[]}));

    //     box_plots_stats(relayoutData,figure,{"stats":[]});

    }
    return dash_clientside.no_update;
}
window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        test_clientside_callback: function(layout,elements) {
            console.log(elements);
            // metanodes = new Map();
            // for(var i =0;i<elements.length;++i){
            //     if("Signatures" in elements[i]["data"]) {
            //         for(var j =0;j<elements[i]["data"]["Signatures"].length;++j){
            //             if(metanodes.has(elements[i]["data"]["Signatures"][j])){
            //                 metanodes.get(elements[i]["data"]["Signatures"][j]).push(elements[i]["position"])
            //             }else{
            //                 metanodes.set(elements[i]["data"]["Signatures"][j],[elements[i]["position"]])
            //             }
            //         }
            //     }
            // }
            // metanodes.forEach((v,k,m) => draw_bounding_box(k,v))
        },
        activate_tooltip: function (mouseoverNodeData,mouseoverEdgeData,fake_graph_size,elements,extent,stylesheets,pos_store_data,cur_children,cur_show,cur_bbox,cur_direction){
            return tooltip(mouseoverNodeData,mouseoverEdgeData,fake_graph_size,elements,extent,stylesheets,pos_store_data,cur_children,cur_show,cur_bbox,cur_direction);
        }
        ,layout_overview:layout_overview
        ,highlight_pathway_neighbourhood:highlight_pathway_neighbourhood
        ,use_elem_weight:use_elem_weight,
        box_plots_stats:box_plots_stats,
        on_box_plot_click:on_box_plot_click,
        display_legend:display_legend
    }
});
function export_image_event_handler(e){
    let graph = e.originalEvent.target.parentNode.parentNode;
    let canvas = graph.querySelector("canvas:nth-of-type(3)");
    switch (graph.id){
        case "overview_graph":
            draw_offscreen(graph.id,"overview.png")
            // download_canvas_image(canvas,"overview.png");
            break;
        case "mono_graph":
            let title = document.getElementById("mono_tab").querySelector("a").innerText;
            draw_offscreen(graph.id,title+".png",document.getElementById("mono_canvas"))
            // download_canvas_image(canvas,title+".png",document.getElementById("mono_canvas"));
            break;
        case "detail_graph":
            draw_offscreen(graph.id,"multi_signature_view.png",document.getElementById("multi_canvas"))
            // download_canvas_image(canvas,"multi_signature_view.png",document.getElementById("multi_canvas"));
            break;
    }
}
function export_text_event_handler(e){
    let graph = e.originalEvent.target.parentNode.parentNode;
    switch (graph.id){
        case "overview_graph":
            export_overview_data();
            break;
        case "mono_graph":
        case "detail_graph":
            dash_clientside.set_props("exportImage", {"value":graph.id});
            document.getElementById("export_json_btn").click();
            break;
    }
}
function node_event_handler(e){
    let href = "";
    if(e.target.data()["is_metanode"]){
        href = e.target.data().tooltip_content.props.children[1].props.href;

    }else{
        href = e.target.data().tooltip_content[1].props.href;
    }
    
    let a = document.createElement("a");
    a.href=href;
    a.target="_blank";
    a.click();
}
window.dashCytoscapeFunctions = Object.assign(
    {},
    window.dashCytoscapeFunctions,
    {
        export_image_event_handler:export_image_event_handler,
        export_text_event_handler:export_text_event_handler,
        node_event_handler:node_event_handler
    })
function ajust_flex(y){
    var upGrow = parseFloat(document.getElementById("overview_col").parentNode.style["flex-grow"]);
    var upBasis = parseFloat(document.getElementById("overview_col").parentNode.style["flex-basis"].split("%"))/100;
    var downBasis = parseFloat(document.getElementById("second_row").style["flex-basis"].split("%"))/100;
    var downGrow = parseFloat(document.getElementById("second_row").style["flex-grow"]);
    var titleHeight = document.getElementById("title").clientHeight;
    var totalHeight = document.querySelector("#full_row > .col-10").clientHeight;
    var upHeight = document.getElementById("overview_col").parentNode.clientHeight;
    var downHeight = document.getElementById("second_row").clientHeight;
    var newUpHeight = y;
    if(Math.abs(newUpHeight-upHeight)>10){
        var newDownHeight = downHeight + upHeight-newUpHeight;
        var freespace = totalHeight-upBasis*totalHeight-downBasis*totalHeight-titleHeight;
        var fgUp = (newUpHeight/*-upBasis*totalHeight-downBasis*totalHeight-titleHeight*/)/(freespace);
        var fgDown = (newDownHeight/*-upBasis*totalHeight-downBasis*totalHeight-titleHeight*/)/(freespace);
        document.getElementById("second_row").style["flex-grow"] = fgDown
        document.getElementById("overview_col").parentNode.style["flex-grow"]=fgUp;
        return true;
    }
    return false;
}

window.addEventListener("resize",function(e){
    var overview_graph = document.getElementById("overview_graph");
    if(overview_graph!==undefined){
        // var elem = document.getElementById('detail_graph');
        // var span = document.getElementById('detail_resize_span');
        // var resize_state = {
        //     "width": {
        //       "is_resizing": false,
        //     },
        //     "height": {
        //       "is_resizing": false
        //     },
        //     "window_resize":true

        //   }
        // dash_clientside.set_props("resize_state",resize_state);
        // dash_clientside.set_props("fake_graph_size",{'width':elem.clientWidth,'height':elem.clientHeight,'width_span':span.clientWidth,"AR":elem.clientWidth/elem.clientHeight});
        // console.log("relayout",{'width':elem.clientWidth,'height':elem.clientHeight,'width_span':span.clientWidth,"AR":elem.clientWidth/elem.clientHeight});
        
        // overview_graph['_cyreg']["cy"].layout({"name":"cose","nodeDimensionsIncludeLabels":true,"animate":false}).run();
        var e2 = new MouseEvent("click");
        document.getElementById("overview_col").dispatchEvent(e2);
}
})
function scale_box_plot_for_stats(h,nb){
    const plot = document.querySelector("#activation_boxplot div.js-plotly-plot");
    let y_max = get_y_max(plot);

    let new_ymax = y_max*(1+h*(nb+1));
    if(new_ymax>plot.layout.yaxis.range[1]){
        Plotly.relayout(plot,{"yaxis.range[1]":new_ymax});
    }


}
function createStatLink(traceIdI,boxIdI,traceIdJ,boxIdJ,h,k,label,shapes){
    const plot = document.querySelector("#activation_boxplot div.js-plotly-plot");
    let { x0, x1 } = stat_shape_x(traceIdI, plot, boxIdI, traceIdJ, boxIdJ);
    if(x0!=undefined && x1 !=undefined){
        let y_max = get_y_max(plot);
        let y = create_shape(y_max, h, k, x0, x1, label, shapes);
        return y+y_max*h;
    }
    return 0;
}
function download_canvas_image(canvas,filename,legend_canvas){
    let oc = new OffscreenCanvas(canvas.width,canvas.height);
    let inCtx = canvas.getContext("2d");
    
    let outCtx = oc.getContext("2d");
    if(legend_canvas!=undefined){
        let inLegendCtx = legend_canvas.getContext("2d");
        outCtx.putImageData(inLegendCtx.getImageData(0,0,legend_canvas.width,legend_canvas.height),0,0);
        createImageBitmap(inCtx.getImageData(0,0,canvas.width,canvas.height)).then(imageBitMap => {
            outCtx.drawImage(imageBitMap,0,0);
            oc.convertToBlob({"type":"image/png"}).then((blob) => {
                let img =URL.createObjectURL(blob);
                //oc.toDataURL("image/png");
                let a = document.createElement("a");
                a.download= filename;
                a.href=img;
                a.click();
            });
        })
        
        
    }else{
        outCtx.putImageData(inCtx.getImageData(0,0,canvas.width,canvas.height),0,0);
        oc.convertToBlob({"type":"image/png"}).then((blob) => {
            let img =URL.createObjectURL(blob);
            //oc.toDataURL("image/png");
            let a = document.createElement("a");
            a.download= filename;
            a.href=img;
            a.click();
        });
    }
    
    
}
function download_plotly_image(plot,filename){
    Plotly.downloadImage(plot,{height:plot.clientHeight,width:plot.clientWidth,format:"png",filename:filename})
}
function stat_shape_x(plot,traceIdI,  boxIdI, traceIdJ, boxIdJ) {
    let parentG = document.querySelector("#activation_boxplot svg.main-svg g.plot g.boxlayer.mlayer");
    let traceIndexI = undefined;
    let boxIndexI = undefined;
    let boxIndexJ = undefined;
    let traceIndexJ = undefined;
    let Gi = undefined;
    let Gj = undefined;

    for (var i = 0; i < parentG.children.length && (traceIdI == undefined || traceIndexJ == undefined); i++) {
        if (parentG.children[i]["__data__"][0].trace.name == traceIdI) {
            traceIndexI = i;
            Gi = parentG.children[traceIndexI];
            for (var j = 0; j < Gi.children.length && (boxIndexI == undefined); j++) {
                let x = Gi.children[j]["__data__"].x;
                if (plot._fullLayout.xaxis.c2d(x) == boxIdI) {
                    boxIndexI = j;
                }
            }

        }
        if (parentG.children[i]["__data__"][0].trace.name == traceIdJ) {
            traceIndexJ = i;
            Gj = parentG.children[traceIndexJ];
            for (var j = 0; j < Gj.children.length && (boxIndexJ == undefined); j++) {
                let x = Gj.children[j]["__data__"].x;
                if (plot._fullLayout.xaxis.c2d(x) == boxIdJ) {
                    boxIndexJ = j;
                }
            }
        }
    }
    let x0 = undefined;
    let x1 = undefined;
    if (traceIndexI != undefined && traceIndexJ != undefined) {
        let boxi = Gi.children[boxIndexI];
        let boxj = Gj.children[boxIndexJ];
        x0 = boxi["__data__"].x + boxi["__data__"].t.bPos;
        x1 = boxj["__data__"].x + boxj["__data__"].t.bPos;
    }
    return x0<x1?{x0:x0,x1:x1 }:{ x0:x1,x1:x0  };
}

function create_shape(y_max, h, k, x0, x1, label, shapes) {
    let y = y_max * (1 + h * (k + 1));
    let line = { type: 'line', x0: x0, y0: y, x1: x1, y1: y, line: { color: 'rgb(0,0,0)', width: 1 }, label: { text: label, font: { color: "black", size: 10 }, textposition: "middle", yanchor: "middle" } };
    let path = { type: 'path', path: `M${x0},${y - h * y_max / 2}V${y}H${x1}V${y - h * y_max / 2}`, line: { color: 'rgb(0,0,0)', width: 1 }, label: { text: label, font: { color: "black" }, textposition: "top center", yanchor: "middle" } };
    shapes.push(path);
    return y;
}

function get_y_max(plot) {
    let y_max = parseFloat("nan");
    const calcdata = plot.calcdata;
    for (var traceId = 0; traceId < calcdata.length; traceId += 1) {
        for (var boxId = 0; boxId < calcdata[traceId].length; boxId += 1) {
            if (calcdata[traceId][boxId].max!==undefined && !(y_max > calcdata[traceId][boxId].max)) {

                y_max = calcdata[traceId][boxId].max;
            }
        }
    }
    return y_max;
}

