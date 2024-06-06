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

function tooltip(mouseoverNodeData,fake_graph_size,elements,extent,stylesheets,pos_store_data,cur_children,cur_show,cur_bbox,cur_direction){
    console.log(mouseoverNodeData);
    if(mouseoverNodeData==null){
    
    //TODO keep cur if close else reset
        //return [cur_children,cur_show,cur_bbox,cur_direction];
        return [[],false,{},"right"];
    }else{
        const cy_elem = document.getElementById(dash_clientside.callback_context.triggered_id)['_cyreg']["cy"].elements("#" + mouseoverNodeData["id"]).toArray();

        var direction="left";
        var elem = undefined;
        console.log(cy_elem)
        if(cy_elem[0].isEdge()){
            console.log(cy_elem[0].isEdge(),cy_elem[0].renderedMidpoint(),cy_elem[0].midpoint(),cy_elem[0].data("symbols"))
            content=cy_elem[0].data("symbols")
            rmp = cy_elem[0].renderedMidpoint()
            mp = cy_elem[0].midpoint()
            
            return [content,true,{"x0":rmp['x'],"y0":rmp['y'],"x1":rmp['x'],"y1":rmp['y']},direction];
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
            var detail_graph = document.getElementById('detail_graph');
            var detail_resize_span = document.getElementById('detail_resize_span');
            //console.log("update_width_start",dash_clientside.callback_context.triggered_id,Object.assign({},fake_graph_size,{'width':elem.clientWidth,'height':elem.clientHeight,'width_span':detail_resize_span.clientWidth,"AR":elem.clientWidth/elem.clientHeight}))
            var x = 0;
            var y=0;
            var width=0;
            var height=0;
            var wratio=detail_graph.clientWidth/(extent["x2"]-extent["x1"]);
            var hratio=detail_graph.clientHeight/(extent["y2"]-extent["y1"]);
            // x = (elem["position"]["x"]-extent["x1"]) *wratio+detail_resize_span.clientWidth;
            // y = (elem["position"]["y"]-extent["y1"]) *hratio
            x = cy_elem[0].renderedPosition("x")+document.getElementById("detail_resize_span").clientWidth
            y = cy_elem[0].renderedPosition("y")
            console.log(elem,cy_elem[0].renderedPosition(),cy_elem[0].relativePosition(),x,y);

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
                if(y_crosses.length<2){
                    for(var i =0; i<polygon.length;i++){
                        if(polygon[i][0]!=0 && polygon[(i+1)%polygon.length][0]!=0){
                            if(polygon[i][1]*polygon[(i+1)%polygon.length][1]<0){
                                y_crosses.push(
                                    polygon[i][0]+(Math.abs(polygon[i][1])/(Math.abs(polygon[i][1])+Math.abs(polygon[(i+1)%polygon.length][1])))*(polygon[(i+1)%polygon.length][0]-polygon[i][0])
                                );
                            }
                        }
                    }
                }
                if(y_crosses.length==2){
                    if(direction=="right"){
                        x+=0.5*y_crosses.reduce((a,b) => Math.max(a,b),-1)*width*wratio;
                    }else{
                        x+=0.5*y_crosses.reduce((a,b) => Math.min(a,b),1)*width*wratio;
                    }
                    width=0;
                    height=0;
                }else{
                    console.log("y_crosses",y_crosses);
                }
            }
            return [mouseoverNodeData["tooltip_content"],true,{"x0":x-width*wratio/2,"y0":y-height*hratio/2,"x1":x+width*wratio/2,"y1":y+height*hratio/2},direction];
        }
    }
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
        activate_tooltip: function (mouseoverNodeData,fake_graph_size,elements,extent,stylesheets,pos_store_data,cur_children,cur_show,cur_bbox,cur_direction){
            return tooltip(mouseoverNodeData,fake_graph_size,elements,extent,stylesheets,pos_store_data,cur_children,cur_show,cur_bbox,cur_direction);
        }
    }
});


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