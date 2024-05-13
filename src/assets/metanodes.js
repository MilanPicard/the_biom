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