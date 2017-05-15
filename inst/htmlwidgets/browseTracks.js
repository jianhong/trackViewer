HTMLWidgets.widget({

  name: 'browseTracks',

  type: 'output',

  factory: function(el, width, height) {
    var svg = d3.select(el).append("svg");
    var menu = d3.select(el).append("input")
                 .attr("type", "button")
                 .attr("name", "export")
                 .attr("id", "export")
                 .attr("value", "exportSVG");
    return {
      renderValue: function(x) {
        //console.log(x);
        function writeDownloadLink(){
            function fireEvent(obj,evt){
              var fireOnThis = obj;
              var evObj;
              if( document.createEvent ) {
                evObj = document.createEvent('MouseEvents');
                evObj.initEvent( evt, true, false );
                fireOnThis.dispatchEvent( evObj );
              } else if( document.createEventObject ) {
                evObj = document.createEventObject();
                fireOnThis.fireEvent( 'on' + evt, evObj );
              }
            }
            if(typeof(ruler)!="undefined"){
                ruler.remove();
                ruler = undefined;
            }
            if(typeof(marginBox)!="undefined"){
                marginBox.remove();
                marginBox = undefined;
            }
            svgAsDataUri(svg.node(), 'trackViewer.svg', function(uri){
                var a = document.createElement('a');
                a.href = uri;
                a.download = 'trackViewer.svg';
                fireEvent(a, 'click');
            });
            draw();
            ruler = new Ruler();
        }
        d3.select("#export")
          .on("click", writeDownloadLink);
          
        svg.attr("width", width)
           .attr("height", height);
        var margin = {top: 20, right: 20, bottom: 30, left: 80};
        x.opacity = [];
        x.fontsize = [];
        for(var k=0; k<x.name.length; k++){
            x.opacity[x.name[k]] = 1;
            x.fontsize[x.name[k]] = x.type[x.name[k]] === "gene" ? 12 : 10;
        }
        var widthF = function() {
            return(+svg.attr("width") - margin.left - margin.right);
        };
        var heightF = function() {
            return(+svg.attr("height") - margin.top - margin.bottom);
        };
        var g;
        var resizeBtn;
        var marginBox;
        var xscale = d3.scaleLinear().domain([x.start, x.end]).rangeRound([0, widthF()]);
        var trackNames = function(){
            return(x.name);
        }
        var redraw = false;
        var editable_ele;
        var parent;
        var bg;
        var color = "#000";
        var cp; //colorPicker;
        var currentId;
        var coords;
        var coor;
        var SVGRect;
        var resize_bg;
        var xHeight=function(){
            var xH=[];
            xH[0] = 0;
            for(var k=0; k<x.name.length; k++){
                xH[k+1] = x.height[trackNames()[k]]*heightF() + xH[k];
            }
            return(xH);
        };
        
        var clearFun = function(){
            if(typeof(editable_ele)!="undefined"){
                if(typeof(bg)!="undefined"){
                    bg.remove();
                    bg = undefined;
                }
                if(typeof(resize_bg)!="undefined"){
                    resize_bg.remove();
                    resize_bg = undefined;
                }
                d3.select("body").on("keypress", null);
                editable_ele.attr("fill", color);
                editable_ele = undefined;
            }
        };
        
        var changeTrackName = function(k, txt){
            var xkeys = d3.keys(x);
            var known = ["name", "chromosome", "start", "end", "strand"];
            xkeys = xkeys.filter(function(el){
                return known.indexOf(el)<0;
            });
            for(var i=0; i<xkeys.length; i++){
                x[xkeys[i]][txt] = x[xkeys[i]][x.name[k]];
                x[xkeys[i]][x.name[k]] = 'undefined';
            }
            x.name[k] = txt;
            
        };
        var changeText = function(txt){
            if(typeof(editable_ele)!="undefined" && typeof(currentId)!="undefined"){
                editable_ele.text(txt);
                changeTrackName(currentId, txt);
            }
        };
        
        var addBg = function(){
            if(typeof(parent)!="undefined" && 
               typeof(SVGRect)!="undefined" && 
               typeof(editable_ele)!="undefined"){
                bg = parent.insert("rect", ":first-child")
                        .attr("x", SVGRect.x)
                        .attr("y", SVGRect.y)
                        .attr("width", SVGRect.width)
                        .attr("height", SVGRect.height)
                        .attr("fill", "#B3D8FD")
                        .attr("transform", editable_ele.attr("transform"));
                editable_ele.attr("fill", "black");
            }
        };
        
        var editText = function(){
            if(typeof(editable_ele)!="undefined"){
            addBg();
            make_resizable();
            var flag = 1;
            d3.select("body").on("keypress", function(){
                // IE fix
                if (!d3.event) d3.event = window.event;
                var e = d3.event;
                //console.log(e);
                switch(e.keyCode){
                    case 13:
                    case 27:
                        clearFun();
                        break;
                    case 8:
                        changeText(editable_ele.text().substring(0, editable_ele.text().length - 1));
                        break;
                    default:
                        if(flag){
                            changeText(e.key);
                            flag = 0;
                        } else {
                            changeText(editable_ele.text()+e.key);
                        }
                }
            });
            }
        };
        
        var removeTarget = function(){
            editable_ele.attr("stroke-width", 2);
            var opacity = editable_ele.style("opacity")===1 ? 0 : 1;
            //console.log(opacity);
            if (confirm("toggle the selection?") === true) {
                editable_ele.style("opacity", opacity);
                if(x.type[trackNames()[currentId]]==="gene"){
                    d3.selectAll("#arrow"+currentId).style("opacity", opacity);
                }
                x.opacity[trackNames()[currentId]] = opacity;
            }
            editable_ele.attr("stroke-width", opacity===0 ? 10 : 1);
        };
        
        var ColorPicker = function (datId=0) {
            if(typeof(cp)!="undefined"){
                cp.remove();
                cp = undefined;
            }
            if(typeof(coor)==="undefined" || 
               typeof(color)==="undefined"  ||
               typeof(editable_ele)==="undefined" ||
               typeof(currentId)==="undefined"){
                return(null);
            }
            var self = this;
            var colorScale = ["#FFD300", "#FFFF00", "#A2F300", "#00DB00", "#00B7FF", 
                              "#1449C4", "#4117C7", "#820AC3", "#DB007C", "#FF0000", 
                              "#FF7400", "#FFAA00"];
            var getColor = function (i) {
                return colorScale[i];
            };
            defaultColor = color || getColor(0);
        
            self.pickedColor = defaultColor;
            self.picked = function (col) {
                //d3.selectAll(".track"+currentId).attr("fill", col);
                if(x.type[trackNames()[currentId]]==="data"){
                    x.tracklist[trackNames()[currentId]].style.color[datId] = col;
                }else{
                    x.tracklist[trackNames()[currentId]].style.color = col;
                }
                color = col;
                draw();
            };
            var clicked = function () {
                self.picked(self.pickedColor);
            };
        
            var pie = d3.pie().sort(null);
            var arc = d3.arc().innerRadius(25).outerRadius(50);
            coor[0] = coor[0]+50-margin.left;
            coor[1] = coor[1]+50-margin.top;
            if(coor[0] > widthF() - 50){
                coor[0] = widthF() - 50;
            }
            if(coor[1] > heightF() - 50){
                coor[1] = heightF() - 50;
            }
            cp = g
                .append("g")
                .attr("width", 100)
                .attr("height", 100)
                .attr("transform", "translate(" + coor[0] +  " " + coor[1] + ")");
            var defaultPlate = cp.append("circle")
                    .attr("fill", defaultColor)
                    .attr("stroke", "#fff")
                    .attr("stroke-width", 2)
                    .attr("r", 10)
                    .attr("cx", 45)
                    .attr("cy", 45)
                    .on("mouseover", function () {
                        var fill = d3.select(this).attr("fill");
                        self.pickedColor = fill;
                        plate.attr("fill", fill);
                    })
                    .on("click", clicked);
            var closePlate = cp.append("g")
                    .attr("width", 20)
                    .attr("height", 20)
                    .attr("transform", "translate(45 -45)");
            closePlate.append("circle")
                    .attr("fill", "#fff")
                    .attr("stroke", "#000")
                    .attr("stroke-width", 1)
                    .attr("r", 10)
                    .attr("cx", 0)
                    .attr("cy", 0)
                    .on("click", function(){
                        cp.remove();
                    });
            closePlate.append("text")
                    .attr("fill", "#000")
                    .attr("x", -5)
                    .attr("y", 5)
                    .text("X")
                    .style("cursor", "default")
                    .on("click", function(){
                        cp.remove();
                    });
        
            var plate = cp.append("circle")
                .attr("fill", defaultColor)
                .attr("stroke", "#fff")
                .attr("stroke-width", 2)
                .attr("r", 25)
                .attr("cx", 0)
                .attr("cy", 0)
                .on("click", clicked);
        
            cp.datum([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
                .selectAll("path")
                .data(pie)
                .enter()
                .append("path")
                .attr("fill", function (d, i) {
                    return getColor(i);
                })
                .attr("stroke", "#fff")
                .attr("stroke-width", 2)
                .attr("d", arc)
                .on("mouseover", function () {
                    var fill = d3.select(this).attr("fill");
                    self.pickedColor = fill;
                    plate.attr("fill", fill);
                })
                .on("click", clicked);
        };
        
        
        function dragstarted(d) {
          d3.select(this).style("cursor", "move").raise().classed("active", true);
          currentId = Number(d3.select(this).attr("kvalue"));
          clearFun();
        }
        function dragended(d) {
          d3.select(this).style("cursor", "default").classed("active", false);
          if(redraw){
              draw();
              redraw = false;
          }
        }
        
        function draggedText(d) {
          var tx = d3.select(this);
          if(/rotate/.exec(tx.attr("transform"))){
              tx.attr("y", d3.event.x).attr("x", -d3.event.y);
          }else{
              tx.attr("x", d3.event.x).attr("y", d3.event.y);
          }
        }
        function draggedCircle(d) {
          d3.select(this).attr("cx", d.x = d3.event.x).attr("cy", d.y = d3.event.y);
        }
        var checkTrackHeight = function(y){
            var xH = xHeight();
            var swap = currentId;
            if(y>xH[xH.length-1]){
                swap = xH.length - 2;
            }else{
                for(var i=0; i<xH.length; i++){
                    if(y<xH[i]){
                        swap = i-1;
                        break;
                    }
                }
            }
            if(swap!=currentId){
                var swapTrack = d3.select("#track"+swap);
                var currentTrack = d3.select("#track"+currentId);
                var ts = xH[currentId];
                swapTrack.attr("transform", "translate(0," + ts + ")");
                //swap x.name
                var tmp = x.name[currentId];
                x.name[currentId] = x.name[swap];
                x.name[swap] = tmp;
                //swap id
                swapTrack.attr("id", "track"+currentId);
                currentTrack.attr("id", "track"+swap);
                currentId = swap;
            }
            redraw = true;
        };
        function draggedTrack(d){
            var tk = d3.select(this);
            tk.attr("transform", "translate(0," + d3.event.y + ")");
            checkTrackHeight(d3.event.y);
        }
        
        var make_resizable = function(){
            if(typeof(editable_ele)==="undefined" ||
               typeof(currentId)==="undefined" ||
               typeof(SVGRect)==="undefined" ||
               typeof(bg)==="undefined" ||
               typeof(parent)==="undefined"){
                return(null);
            }
            
            var handlerRadius = 2.5;
            resize_bg = parent.append("g");
            var gHandlerPoints = resize_bg.append("circle")
                          .attr("class", "handlerPoint")
                          .style("fill", "white")
                          .style("stroke", "blue")
                          .attr("stroke","")
                          .attr("transform", "translate("+(SVGRect.x + SVGRect.width/2)+","+SVGRect.y+")")
                          .attr("r", handlerRadius)
                          .call(d3.drag()
                                      .on("start", resize_dragstarted)
                                      .on("drag", resize_draged)
                                      .on("end", resize_dragended));
            if(/rotate/.exec(editable_ele.attr("transform"))){
                gHandlerPoints.attr("transform", "rotate(-90) translate("+(SVGRect.x + SVGRect.width)+","+(SVGRect.y+SVGRect.height/2)+")");
            }
            function resize_draged(d){
                var fontsize = x.fontsize[x.name[currentId]] - d3.event.dy;
                editable_ele.style("font-size", fontsize + "px");
                x.fontsize[x.name[currentId]] = fontsize;
                SVGRect = editable_ele.node().getBBox();
                gHandlerPoints.attr("transform", 
                    /rotate/.exec(editable_ele.attr("transform")) ? 
                        "rotate(-90) translate("+(SVGRect.x + SVGRect.width)+","+(SVGRect.y+SVGRect.height/2)+")" : 
                        "translate("+(SVGRect.x + SVGRect.width/2)+","+SVGRect.y+")");
                bg.attr("x", SVGRect.x)
                  .attr("y", SVGRect.y)
                  .attr("width", SVGRect.width)
                  .attr("height", SVGRect.height);
            }
            function resize_dragstarted(d) {
              d3.select(this).style("cursor", "ns-resize").raise().classed("active", true);
            }
        
            function resize_dragended(d) {
              d3.select(this).style("cursor", "default").classed("active", false);
              resize_bg.remove();
              resize_bg = undefined;
              draw();
            }
        };
        
        var make_editable = function(){
            //remove other editalbe
            clearFun();
            // set current editable element info;
            editable_ele = d3.select(this);
            parent = d3.select(this.parentNode);
            currentId = Number(editable_ele.attr("kvalue"));
            coords = d3.mouse(this);
            coor = [d3.event.pageX, d3.event.pageY];
            color = editable_ele.attr("fill");
            SVGRect = this.getBBox();
            var cls = editable_ele.attr("class");
            if(/geneSymbol/.exec(cls)){
                //console.log("symbol");
                editText();
            }
            if(/dataYlabel/.exec(cls)){
                //console.log("dataYlabel");
                editText();
            }
            if(/dataPath1/.exec(cls)){
                //console.log("dataPath1");
                ColorPicker(0);
            }
            if(/dataPath2/.exec(cls)){
                //console.log("dataPath2");
                ColorPicker(1);
            }
            if(/exon/.exec(cls)){
                //console.log("exon");
                ColorPicker(0);
            }
            if(/dataBaseline/.exec(cls)){
                //console.log("dataBaseline");
                removeTarget();
            }
            if(/geneArrow/.exec(cls)){
                //console.log("geneArrow");
                removeTarget();
            }
            if(/geneBaseline/.exec(cls)){
                //console.log("geneBaseline");
                removeTarget();
            }
        };
        
        var dataTrack = function(layer, track, start, end, xscale, yscale, line, k){
            var color = track.style.color;
            //console.log(color);
            //signal dat
            if(track.dat.length>0){
                var data=[{"x":start-1, "y":0}];
                for(var i=0; i<track.dat.length; i++){
                    data.push({
                        "x" : i+start, 
                        "y" : track.dat[i]
                    });
                }
                data.push({"x":end+1, "y":0});
                //console.log(data);
                layer.append("path")
                    .datum(data)
                    .attr("fill", color[0])
                    .attr("stroke", "none")
                    .attr("d", line)
                    .attr("class", "dataPath1 track"+k)
                    .attr("kvalue", k)
                    .on("click", make_editable);
            }
            //signal dat2
            if(track.dat2.length>0){
                var data=[{"x":start-1, "y":0}];
                for(var i=0; i<track.dat2.length; i++){
                    data.push({
                        "x" : i+start, 
                        "y" : -track.dat2[i]
                    });
                }
                data.push({"x":end+1, "y":0});
                //console.log(data);
                layer.append("path")
                    .datum(data)
                    .attr("fill", color[1])
                    .attr("stroke", "none")
                    .attr("d", line)
                    .attr("class", "dataPath2 track"+k)
                    .attr("kvalue", k)
                    .on("click", make_editable);
            }
            // add a line at 0
            layer.append("line")
                 .attr("x1", xscale(start))
                 .attr("x2", xscale(end))
                 .attr("y1", yscale(0))
                 .attr("y2", yscale(0))
                 .attr("stroke", color[0])
                 .attr("stroke-width", x.opacity[x.name[k]]===1 ? "1px" : "10px")
                 .attr("opacity", x.opacity[x.name[k]])
                 .attr("class", "dataBaseline track"+k)
                 .attr("kvalue", k)
                 .on("click", make_editable);
        };
        
        var plotArrow = function(container, data, start, end, yscale, strand, color, k){
            if(strand==="*") return(null);
            var marker = container.append("marker")
                     .attr("id", "arrow"+k)
                     .attr("kvalue", k)
                     .attr("viewBox", "0, 0, 10, 10")
                     .attr("refX", 0)
                     .attr("refY", 5)
                     .attr("markerUnits", "strokeWidth")
                     .attr("markerWidth", 8)
                     .attr("markerHeight", 8)
                     .attr("fill", "none")
                     .attr("stroke", color)
                     .attr("opacity", x.opacity[x.name[k]])
                     .attr("orient", strand==="-"?0:180)
                     .append("path")
                     .attr("d", "M 0 0 L 10 5 L 0 10");
            for(var i=0; i<data.length; i++){
                for(var j=data[i].x0; j<=data[i].x1; j+=10){
                    if(j > start && j < end){
                        container.append("line")
                            .attr("x1",  j)
                            .attr("y1", yscale(0.5))
                            .attr("x2",  j)
                            .attr("y2", yscale(0.5))
                            .attr("marker-end", "url(#arrow"+k+")")
                            .attr("stroke-width", 1)
                            .attr("class", "geneArrow track"+k)
                            .attr("kvalue", k)
                            .on("click", make_editable);
                    }
                }
            }
            
        }
        
        var geneTrack = function(layer, track, start, end, xscale, yscale, wscale, line, k, label){
            //console.log(track);
            var color = track.style.color;
            if(track.dat.start.length>0){
                var data=[];
                var intron=[];
                var geneStart=end, geneEnd=start;
                var intronStart=end+5, intronEnd=start-5;
                var Y=1;
                var feature={
                    "utr5" : Y/2,
                    "CDS"  : Y,
                    "utr3" : Y/2
                };
                for(var i=0; i<track.dat.start.length; i++){
                    if(i > 0){
                        intronStart = track.dat.end[i-1] + 1;
                        intronEnd = track.dat.start[i] - 1;
                        if(intronStart < intronEnd){
                            intron.push({
                                "x0" : xscale(intronStart) + 5,
                                "x1" : xscale(intronEnd)
                            });
                        }
                    }
                    var clipx0 = track.dat.start[i];
                    var clipx1 = track.dat.end[i];
                    if(clipx0<=end && clipx1>=start){
                        if(clipx0<start){
                            clipx0 = start;
                        }
                        if(clipx1>end){
                            clipx1 = end;
                        }
                        data.push({
                            "x" : clipx0,
                            "h" : feature[track.dat.feature[i]],
                            "w" : clipx1 - clipx0 + 1
                        });
                    }
                    if(track.dat.start[i]<geneStart){
                        geneStart = track.dat.start[i];
                    }
                    if(track.dat.end[i]>geneEnd){
                        geneEnd = track.dat.end[i];
                    }
                }
                // add a center line
                if(geneStart < start) geneStart = start;
                if(geneEnd > end) geneEnd = end;
                if(geneStart < geneEnd){
                    layer.append("line")
                         .attr("x1", xscale(geneStart))
                         .attr("x2", xscale(geneEnd))
                         .attr("y1", yscale(0.5))
                         .attr("y2", yscale(0.5))
                         .attr("stroke", color)
                         .attr("stroke-width", x.opacity[label]===1 ? "1px" : "10px")
                         .attr("opacity", x.opacity[label])
                         .attr("class", "geneBaseline track"+k)
                         .attr("kvalue", k)
                         .on("click", make_editable);
                    // add arrows to center line
                    plotArrow(layer, intron, xscale(start), xscale(end), yscale, track.strand, color, k);
                }
                //console.log(intron);
                layer.selectAll(".exon"+k)
                    .data(data).enter().append("rect")
                    .attr("class", "exon track"+k)
                    .attr("kvalue", k)
                    .attr("x", d=> xscale(d.x))
                    .attr("y", d=> yscale(0.5 + d.h/2))
                    .attr("width", d=> wscale(d.w))
                    .attr("height", d=> yscale(1 - d.h))
                    .attr("fill", color)
                    .on("click", make_editable);
                // add gene symbols
                layer.append("g").append("text")
                     .attr("fill", color)
                     .attr("x", xscale(geneStart) - 10)
                     .attr("y", yscale(0.5)+4)
                     .attr("text-anchor", "end")
                     .attr("class", "geneSymbol track"+k)
                     .attr("kvalue", k)
                     .style("font-size", x.fontsize[label]+"px")
                     .text(label)
                     .on("click", make_editable)
                     .call(d3.drag().on("start", dragstarted)
                                    .on("drag", draggedText)
                                    .on("end", dragended));
            }
        };
        
        var xaxis = function(layer, xscale, height){
            layer.append("g")
                    .attr("transform", "translate(0," + height + ")")
                    .call(d3.axisBottom(xscale));
        };
        var yaxis = function(layer, yscale, ylim, height, label, k){
            var ax = layer.append("g")
                .call(d3.axisLeft(yscale).tickSize(5).tickValues(ylim));
            var ylabel = ax
                .append("text")
                .attr("fill", "#000")
                .attr("transform", "rotate(-90)")
                .attr("y", -8)
                .attr("x", -height/2)
                .attr("text-anchor", "middle")
                .attr("class", "dataYlabel track"+k)
                .attr("font-size", x.fontsize[label] + "px")
                .attr("kvalue", k)
                .text(label)
                .on("click", make_editable)
                .call(d3.drag().on("start", dragstarted)
                                    .on("drag", draggedText)
                                    .on("end", dragended));
        };
        var Margin = function(){
            if(typeof(margin)==="undefined"){
                return(null);
            }
            var self = this;
            self.margintop = svg.append("line")
                                 .attr("stroke", "white")
                                 .attr("stroke-width", "3px")
                                 .attr("x1", margin.left)
                                 .attr("y1", margin.top)
                                 .attr("x2", +svg.attr("width")-margin.right)
                                 .attr("y2", margin.top)
                                 .attr("ref", 3)
                                 .style("opacity", 0)
                                 .style("cursor", "ns-resize")
                                 .call(d3.drag().on("drag", draggedMargin));
            self.marginbottom = svg.append("line")
                                 .attr("stroke", "white")
                                 .attr("stroke-width", "3px")
                                 .attr("x1", margin.left)
                                 .attr("y1", +svg.attr("height")-margin.bottom)
                                 .attr("x2", +svg.attr("width")-margin.right)
                                 .attr("y2", +svg.attr("height")-margin.bottom)
                                 .attr("ref", 1)
                                 .style("opacity", 0)
                                 .style("cursor", "ns-resize")
                                 .call(d3.drag().on("drag", draggedMargin));
            self.marginleft = svg.append("line")
                                 .attr("stroke", "white")
                                 .attr("stroke-width", "3px")
                                 .attr("x1", margin.left)
                                 .attr("y1", margin.top)
                                 .attr("x2", margin.left)
                                 .attr("y2", +svg.attr("height")-margin.bottom)
                                 .attr("ref", 2)
                                 .style("opacity", 0)
                                 .style("cursor", "ew-resize")
                                 .call(d3.drag().on("drag", draggedMargin));
            self.marginright = svg.append("line")
                                 .attr("stroke", "white")
                                 .attr("stroke-width", "3px")
                                 .attr("x1", +svg.attr("width")-margin.right)
                                 .attr("y1", margin.top)
                                 .attr("x2",+svg.attr("width")- margin.right)
                                 .attr("y2", +svg.attr("height")-margin.bottom)
                                 .attr("ref", 4)
                                 .style("opacity", 0)
                                 .style("cursor", "ew-resize")
                                 .call(d3.drag().on("drag", draggedMargin));
            self.remove = function(){
                self.margintop.remove();
                self.marginleft.remove();
                self.marginbottom.remove();
                self.marginright.remove();
            }
            
            function draggedMargin(d){
                var coordinates = d3.mouse(svg.node());
                var dy = coordinates[1];
                var dx = coordinates[0];
                switch(Number(d3.select(this).attr("ref"))){
                    case 1:
                        dy = +svg.attr("height") - dy;
                        if(dy<0) dy=0;
                        margin.bottom = dy;
                        self.marginbottom.attr("y1", dy).attr("y2", dy);
                        break;
                    case 2:
                        if(dx<0) dx=0;
                        margin.left = dx;
                        self.marginleft.attr("x1", dx).attr("x2", dx);
                        break;
                    case 3:
                        if(dy<0) dy=0;
                        margin.top = dy;
                        self.margintop.attr("y1", dy).attr("y2", dy);
                        break;
                    case 4:
                        dx = +svg.attr("width") - dx;
                        if(dx<0) dx=0;
                        margin.right = dx;
                        self.marginright.attr("x1", dx).attr("x2", dx);
                        break;
                }
                draw();
            }
            return(self);
        };
        var draw = function(vspace=10){
            if(typeof(g)!="undefined") g.remove();
            if(typeof(resizeBtn)!="undefined") resizeBtn.remove();
            if(typeof(marginBox)!="undefined") marginBox.remove();
            xscale = d3.scaleLinear().domain([x.start, x.end]).rangeRound([0, widthF()]);
            g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
            g.append("rect")
             .attr("x", -margin.left)
             .attr("y", -margin.right)
             .attr("width", svg.attr("width"))
             .attr("height", svg.attr("height"))
             .attr("fill", "white")
             .attr("stroke", "none")
             .style("opacity", 0)
             .on("click", function(){
                clearFun();
            });
            var currH = 0;
            var wscale = d3.scaleLinear().domain([0, x.end-x.start+1]).rangeRound([0, widthF()]);
            
            for(var k=0; k<trackNames().length; k++){
                var thisHeight = x.height[trackNames()[k]]*heightF() - 2*vspace;
                var track = g.append("g")
                              .attr("kvalue", k)
                              .attr("id", "track"+k)
                              .attr("width", widthF())
                              .attr("height", thisHeight)
                              .attr("transform", "translate(0," + (currH + vspace) +")");
                //make track drag-able
                track.call(d3.drag().on("start", dragstarted)
                                    .on("drag", draggedTrack)
                                    .on("end", dragended));
                var currTrack = x.tracklist[trackNames()[k]];
                var yscale = d3.scaleLinear().domain(currTrack.ylim).rangeRound([thisHeight, 0]);
                //line x y;
                var line = d3.line()
                    .x(function(d) {return xscale(d.x);})
                    .y(function(d) {return yscale(d.y);});
                if(x.type[trackNames()[k]] === "data"){
                    //console.log(currTrack);
                    //console.log(currTrack.ylim);
                    dataTrack(track, currTrack, x.start, x.end, xscale, yscale, line, k);
                    // x,y-axis
                    //xaxis(track, xscale, x.height[trackNames()[k]]*height);
                    yaxis(track, yscale, currTrack.ylim, thisHeight, trackNames()[k], k);
                }else{
                    geneTrack(track, currTrack, x.start, x.end, xscale, yscale, wscale, line, k, trackNames()[k]);
                }
                currH +=x.height[trackNames()[k]]*heightF();
            }
            xaxis(g, xscale, heightF());
            
            resizeBtn = svg.append("rect").attr("fill", "white")
                               .attr("x", svg.attr("width")-10)
                               .attr("y", svg.attr("height")-10)
                               .attr("width", 10)
                               .attr("height", 10)
                               .style("cursor", "nwse-resize")
                               .call(d3.drag().on("drag", function(d){
                                svg.attr("width", d3.event.x).attr("height", d3.event.y);
                                if(typeof(ruler)!="undefined") ruler.rulemove(d3.mouse(this));
                                draw();
                               }));
            marginBox = new Margin();
        };
        
        var Ruler = function(){
            var self = this;
            self.xrule = svg.insert("line", ":first-child")
                                 .attr("stroke", "#ccc")
                                 .attr("stroke-width", "1px")
                                 .style("stroke-dasharray", ("3, 3"));
            self.xrule_label = svg.insert("text", ":first-child")
                                    .attr("fill", "#333")
                                    .attr("y", 10)
                                    .attr("x", 8)
                                    .attr("text-anchor", "end")
                                    .style("font-size", "0.75em")
                                    .text("");
            self.yrule = svg.insert("line", ":first-child")
                             .attr("stroke", "#ccc")
                             .attr("stroke-width", "1px")
                             .style("stroke-dasharray", ("3, 3"));
            self.yrule_label = svg.insert("text", ":first-child")
                                    .attr("fill", "#333")
                                    .attr("y", 10)
                                    .attr("x", 8)
                                    .attr("text-anchor", "start")
                                    .style("font-size", "0.75em")
                                    .text("");
            self.trackNum = function(xpos, ypos){
                var H = (ypos - margin.top)/heightF();
                for(var i=0; i<trackNames().length; i++){
                    H = H - x.height[trackNames()[i]];
                    if(H<0){
                        if(x.type[trackNames()[i]]==="data"){
                            var a=x.tracklist[trackNames()[i]].dat,
                                b=x.tracklist[trackNames()[i]].dat2;
                            var c="";
                            if(a.length) c = d3.format(".2f")(a[xpos]);
                            if(b.length) c = c + "; "+d3.format(".2f")(b[xpos]);
                            return(c);
                        }else{
                            return(ypos);
                        }
                    }
                }
                return(i);
            }
            self.rulemove = function(coords){
                var xpos = Math.round(xscale.invert(coords[0]-margin.left));
                var ypos = Math.round(d3.event.pageY);
                self.xrule.attr("x1", coords[0])
                     .attr("x2", coords[0])
                     .attr("y1", 0)
                     .attr("y2", svg.attr("height"));
                self.yrule.attr("x1", 0)
                     .attr("x2", svg.attr("width"))
                     .attr("y1", coords[1])
                     .attr("y2", coords[1]);
                self.xrule_label.attr("y", 10)
                           .attr("x", coords[0]-5)
                           .text(xpos);
                self.yrule_label.attr("y", coords[1]-5)
                           .attr("x", 3)
                           .text(self.trackNum(xpos-x.start, ypos));
            };
            self.remove = function(){
                self.xrule.remove();
                self.yrule.remove();
                self.xrule_label.remove();
                self.yrule_label.remove();
                svg.on("mousemove", null);
            }
            svg.on("mousemove", function(){
                self.rulemove(d3.mouse(this));
            });
            return(self);
        }
        
        draw();
        var ruler = new Ruler();
      },

      resize: function(width, height) {
        svg.attr("width", width)
           .attr("height", height);
      },

      svg: svg
    };
  }
});
