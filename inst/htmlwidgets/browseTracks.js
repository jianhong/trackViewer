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
    var menu2 = d3.select(el).append("input")
                 .attr("type", "button")
                 .attr("name", "exportPNG")
                 .attr("id", "exportPNG")
                 .attr("value", "exportPNG");
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
        d3.select("#exportPNG")
          .on("click", function(){
            if(typeof(ruler)!="undefined"){
                ruler.remove();
                ruler = undefined;
            }
            if(typeof(marginBox)!="undefined"){
                marginBox.remove();
                marginBox = undefined;
            }
            var background = svg.insert("rect", ":first-child")
                              .attr("width", svg.attr("width"))
                              .attr("height", svg.attr("height"))
                              .attr("fill", "white")
                              .attr("strock", "none");
            saveSvgAsPng(svg.node(), "trackViewer.png", {scale: 4});
            background.remove();
            draw();
            ruler = new Ruler();
          });
          
        svg.attr("width", width)
           .attr("height", height);
        var margin = {top: 20, right: 20, bottom: 30, left: 100};
        x.opacity = {};
        x.fontsize = {};
        x.markers = {};
        var defaultFontSize = 16;
        var defaultTickFontSize = 16;
        for(var k=0; k<x.name.length; k++){
            x.opacity[x.name[k]] = 1;
            x.fontsize[x.name[k]] = defaultFontSize;
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
        var markerGroup;
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
            if(x.name.indexOf(txt) != -1){
                console.log("used name" + txt);
                return;
            }
            var xkeys = d3.keys(x);
            var known = ["name", "chromosome", "start", "end", "strand", "markers"];
            xkeys = xkeys.filter(function(el){
                return known.indexOf(el)<0;
            });
            for(var i=0; i<xkeys.length; i++){
                x[xkeys[i]][txt] = x[xkeys[i]][x.name[k]];
                delete(x[xkeys[i]][x.name[k]]);
            }
            x.name[k] = txt;
            //console.log(x);
        };
        var changeText = function(txt){
            if(typeof(editable_ele)!="undefined" && typeof(currentId)!="undefined"){
                editable_ele.text(txt);
                if(/Mlabel/.exec(editable_ele.attr("class"))){
                    x.markers[editable_ele.attr("ref")].txt = txt;
                }else{
                    changeTrackName(currentId, txt);
                }
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
                    case 46:
                        if(/Mlabel/.exec(editable_ele.attr("class"))){
                            delete(x.markers[editable_ele.attr("ref")]);
                            editable_ele.remove();
                            clearFun();
                        }
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
            var opacity = editable_ele.style("opacity")==1 ? 0 : 1;
            //console.log(opacity);
            if (confirm("toggle the selection?") == true) {
                editable_ele.style("opacity", opacity);
                if(x.type[trackNames()[currentId]]==="gene"){
                    d3.selectAll("#arrow"+currentId).style("opacity", opacity);
                }
                x.opacity[trackNames()[currentId]] = opacity;
            }
            editable_ele.attr("stroke-width", opacity==0 ? 10 : 1);
        };
        
        var ColorPicker = function (datId=0) {
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
                if(datId==3){
                    //markers
                    x.markers[editable_ele.attr("ref")].color = col;
                }else{
                    if(x.type[trackNames()[currentId]]==="data"){
                        x.tracklist[trackNames()[currentId]].style.color[datId] = col;
                    }else{
                        x.tracklist[trackNames()[currentId]].style.color = col;
                    }
                }
                color = col;
                draw(); // refresh the tracks.
                newCP();// keep it on top
            };
            var clicked = function () {
                self.picked(self.pickedColor);
            };
        
            var pie = d3.pie().sort(null);
            var arc = d3.arc().innerRadius(25).outerRadius(50);
            var currentCoor = coor;
            currentCoor[0] = currentCoor[0]+50;
            currentCoor[1] = currentCoor[1]+50;
            if(currentCoor[0] > widthF() - 50){
                currentCoor[0] = widthF() - 50;
            }
            if(currentCoor[1] > heightF() - 50){
                currentCoor[1] = heightF() - 50;
            }
            var newCP = function(){
                if(typeof(cp)!="undefined"){
                    cp.remove();
                    cp = undefined;
                }
                cp = svg
                .append("g")
                .attr("width", 100)
                .attr("height", 100)
                .attr("transform", "translate(" + currentCoor[0] +  " " + currentCoor[1] + ")")
                .call(d3.drag().on("drag", function(d){//moveable;
                    currentCoor = [currentCoor[0]+d3.event.dx, currentCoor[1]+d3.event.dy];
                    d3.select(this)
                      .style("cursor", "move")
                      .attr("transform", "translate(" + currentCoor[0] +  " " + currentCoor[1] + ")");
                }).on("end", function(d){
                    d3.select(this).style("cursor", "default");
                }));
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
            var blackPlate = cp.append("circle")
                    .attr("fill", "black")
                    .attr("stroke", "#fff")
                    .attr("stroke-width", 2)
                    .attr("r", 10)
                    .attr("cx", -45)
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
            newCP();
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
                var fontsize = 12;
                if(/Mlabel/.exec(editable_ele.attr("class"))){
                    fontsize = x.markers[editable_ele.attr("ref")].fontsize - d3.event.dy;
                    x.markers[editable_ele.attr("ref")].fontsize = fontsize;
                }else{
                    if(parent.attr("class")==="tick"){
                        defaultTickFontSize = defaultTickFontSize - d3.event.dy;
                        fontsize = defaultTickFontSize;
                    }else{
                        fontsize = x.fontsize[x.name[currentId]] - d3.event.dy;
                        x.fontsize[x.name[currentId]] = fontsize;
                    }
                }
                editable_ele.style("font-size", fontsize + "px");
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
            if(parent.attr("class")==="tick"){
                //console.log("tick");
                addBg();
                make_resizable();
            }
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
            if(/Marker/.exec(cls)){
                //console.log("marker");
                ColorPicker(3);
            }
            if(/Mlabel/.exec(cls)){
                //console.log("marker labels");
                editText();
                ColorPicker(3);
            }
            if(/arrowline/.exec(cls)){
                //console.log("arrow");
                ColorPicker(3);
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
                 .attr("stroke-width", x.opacity[x.name[k]]==1 ? "1px" : "10px")
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
                         .attr("stroke-width", x.opacity[label]==1 ? "1px" : "10px")
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
                    .call(d3.axisBottom(xscale).tickSize(5).ticks(5))
                    .attr("font-size", defaultTickFontSize);
        };
        var yaxis = function(layer, yscale, ylim, height, label, k){
            var ax = layer.append("g")
                .call(d3.axisLeft(yscale).tickSize(5).tickValues(ylim))
                .attr("font-size", defaultTickFontSize);
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
        
        var Marker = function(){
            if(typeof(x.markers)==="undefined"){
                x.markers = [];
            }
            var self = this;
            self.svgDefault = true;
            self.remove = function(){
                if(typeof(self.g)!="undefined") self.g.remove();
                svg.on("dblclick", null);
            };
            var makeLabelDraggable = d3.drag()
                       .on("drag", function(d){
                            var curr = d3.select(this);
                            var coords = d3.mouse(svg.node());
                            var posx = Math.round(xscale.invert(coords[0]-margin.left));
                            var posy = Math.round(coords[1]);
                            var m = x.markers[curr.attr("ref")];
                            if(posx!=m.ref[0] && posy!=m.ref[1]){
                                m.ref = [posx, posy];
                                x.markers["text"+posx + "_" + posy] = m;
                                delete(x.markers[curr.attr("ref")]);
                                curr.attr("id", "text"+m.ref[0] + "_" + m.ref[1])
                                   .attr("ref", "text"+m.ref[0] + "_" + m.ref[1])
                                   .attr("x", xscale(m.ref[0]) + margin.left)
                                   .attr("y", m.ref[1])
                                   .style("cursor", "move");
                            }
                            if(typeof(bg)!="undefined"){
                                bg.remove();
                            }
                            if(typeof(resize_bg)!="undefined"){
                                resize_bg.remove();
                            }
                       })
                       .on("end", function(d){
                            d3.select(this).style("cursor", "default");
                       });
            svg.on("dblclick", function(){
                if(self.svgDefault){
                    var coords = d3.mouse(svg.node());
                    var posx = Math.round(xscale.invert(coords[0]-margin.left));
                    var posy = Math.round(coords[1]); // can not keep position
                    var m = {
                        "ref" : [posx, posy],
                        "markertype" : 2,
                        "color" : "black",
                        "opacity" : 1,
                        "fontsize" : 12,
                        "txt" : "label"
                    };
                    x.markers["text"+posx + "_" + posy] = m;
                    self.g.append("text")
                           .attr("id", "text"+m.ref[0] + "_" + m.ref[1])
                           .attr("ref", "text"+m.ref[0] + "_" + m.ref[1])
                           .attr("fill", m.color)
                           .attr("y", m.ref[1])
                           .attr("x", xscale(m.ref[0]) + margin.left)
                           .attr("text-anchor", "start")
                           .style("font-size", m.fontsize+"px")
                           .text(m.txt)
                           .attr("class", "Mlabel")
                           .on("click", make_editable)
                           .call(makeLabelDraggable);
                    var m1 = {
                        "ref" : [posx, posy, Math.round(xscale.invert(coords[0]-margin.left-20)), posy+20],
                        "markertype" : 3,
                        "color" : "#000",
                        "opacity" : 1,
                        "linetype" : "solid",
                        "linewidth" : 1
                    };
                    x.markers["arrow"+m1.ref[0]+"_"+m1.ref[2]+"_"+m1.ref[1]+"_"+m1.ref[3]]=m1;
                    var arrow = new Arrow(self.g, m1);
                }else{
                    self.svgDefault=true;
                }
            });
            self.g = svg.append("g");
            self.draggedLine = function(d){
                var coords = d3.mouse(svg.node())[0];
                d3.select(this).attr("x1", coords).attr("x2", coords);
            }
            self.dragendLine = function(d){
                var coords = d3.mouse(svg.node());
                var posx = Math.round(xscale.invert(coords[0]-margin.left));
                var old = d3.select(this).attr("ref");
                if(old != "line"+posx){
                    x.markers["line"+posx] = x.markers[old];
                    delete(x.markers[old]);
                    x.markers["line"+posx].ref = posx;
                    d3.select(this).attr("ref", "line"+posx);
                    self.redraw();
                }
            }
            self.dragstarted = function(d){
                            coor = d3.mouse(this);
                            var posx = Math.round(xscale.invert(coor[0] - margin.left));
                            var m = {
                                "ref" : posx,
                                "markertype" : 1,
                                "color" : "gray",
                                "opacity" : 0.3,
                                "linetype" : "solid",
                                "linewidth" : 1
                            };
                            self.g.append("rect")
                                           .attr("id", "rect"+posx)
                                           .attr("stroke", "none")
                                           .attr("fill", m.color)
                                           .attr("x", xscale(m.ref) + margin.left)
                                           .attr("width", m.linewidth)
                                           .attr("y", margin.top)
                                           .attr("height", +svg.attr("height") - margin.bottom - margin.top)
                                           .attr("ref", "rect"+m.ref)
                                           .style("opacity", m.opacity)
                                           .attr("class", "Marker");
                };
              self.dragged = function(d){
                                    var coords = d3.mouse(this);
                                    var posx = Math.round(xscale.invert(coor[0] - margin.left));
                                    d3.select("#rect"+posx).attr("width", coords[0] - coor[0]);
                                    };
              self.dragended = function(d){
                                var coords = d3.mouse(this);
                                var posx = Math.round(xscale.invert(coor[0] - margin.left));
                                if(coords[0] - coor[0] > 3){
                                    x.markers["rect"+posx] = {
                                        "ref" : posx,
                                        "markertype" : 1,
                                        "color" : "gray",
                                        "opacity" : 0.3,
                                        "linetype" : "solid",
                                        "linewidth" : xscale.invert(coords[0] - margin.left) - posx
                                    };
                                }else{
                                    x.markers["line"+posx] = {
                                        "ref" : posx,
                                        "markertype" : 0,
                                        "color" : "gray",
                                        "opacity" : 1,
                                        "linetype" : "dashed",
                                        "linewidth" : 1
                                    };
                                }
                                self.redraw();
                            };
            self.resizeRectL = function(d){
                var coords = d3.mouse(this);
                var posx = Math.round(xscale.invert(coords[0] - margin.left));
                var obj = d3.select(this);
                var ref = obj.attr("ref");
                if(posx != x.markers[ref].ref){
                    x.markers[ref].linewidth = x.markers[ref].ref - posx + x.markers[ref].linewidth;
                    obj.attr("ref", "rect" + posx)
                       .attr("id", "rectL" + posx)
                       .attr("x1", xscale(posx) + margin.left)
                       .attr("x2", xscale(posx) + margin.left);
                    d3.select("#rect"+x.markers[ref].ref)
                      .attr("ref", "rect" + posx)
                      .attr("id", "rect" + posx)
                      .attr("x", xscale(posx) + margin.left)
                      .attr("width", xscale(posx + x.markers[ref].linewidth) - xscale(posx));
                    d3.select("#rectR"+x.markers[ref].ref)
                      .attr("ref", "rect" + posx)
                      .attr("id", "rectR" + posx)
                      .attr("x1", xscale(x.markers[ref].linewidth + posx)+margin.left)
                      .attr("x2", xscale(x.markers[ref].linewidth + posx)+margin.left)
                    x.markers[ref].ref = posx;
                    x.markers["rect"+posx] = x.markers[ref];
                    delete(x.markers[ref]);
                }
            };
            self.resizeRectR = function(d){
                var coords = d3.mouse(this);
                var posx = Math.round(xscale.invert(coords[0] - margin.left));
                var obj = d3.select(this);
                var ref = obj.attr("ref");
                if(posx != x.markers[ref].ref + x.markers[ref].linewidth){
                    obj.attr("x1", xscale(posx) + margin.left)
                       .attr("x2", xscale(posx) + margin.left);
                    d3.select("#rect"+x.markers[ref].ref)
                      .attr("width", xscale(posx) - xscale(x.markers[ref].ref));
                    x.markers[ref].linewidth = posx - x.markers[ref].ref;
                }
            };
            var Arrow = function(arrowcontainer, m){
                var arrowline = this;
                arrowline.remove = function(){
                    arrowline.g.remove();
                };
                var x1=m.ref[0], x2=m.ref[2], y1=m.ref[1], y2=m.ref[3];
                arrowline.g = arrowcontainer.append("g")
                                       .attr("class", "Marker")
                                       .attr("ref", "arrow"+x1+"_"+x2+"_"+y1+"_"+y2);
                var marker = arrowline.g.append("marker")
                         .attr("id", "arrowhead"+x1+"_"+x2+"_"+y1+"_"+y2)
                         .attr("viewBox", "0, 0, 10, 10")
                         .attr("refX", 0)
                         .attr("refY", 5)
                         .attr("markerUnits", "strokeWidth")
                         .attr("markerWidth", 8)
                         .attr("markerHeight", 8)
                         .attr("fill", m.color)
                         .attr("stroke", m.color)
                         .attr("opacity", m.opacity)
                         .attr("orient", "auto")
                         .append("path")
                         .attr("d", "M 0 0 L 10 5 L 0 10 Z");
                arrowline.line = arrowline.g.append("line")
                                    .attr("x1", xscale(x1)+margin.left)
                                    .attr("x2", xscale(x2)+margin.left)
                                    .attr("y1", y1)
                                    .attr("y2", y2)
                                    .attr("stroke", m.color)
                                    .attr("stroke-width", m.linewidth)
                                    .attr("opacity", m.opacity)
                                    .attr("class", "arrowline")
                                    .attr("ref", "arrow"+x1+"_"+x2+"_"+y1+"_"+y2)
                                    .attr("marker-end", "url(#arrowhead"+x1+"_"+x2+"_"+y1+"_"+y2+")");
                arrowline.linecover = arrowline.g.append("line")
                                    .attr("x1", xscale(x1)+margin.left)
                                    .attr("x2", xscale(x2)+margin.left)
                                    .attr("y1", y1)
                                    .attr("y2", y2)
                                    .attr("stroke", m.color)
                                    .attr("stroke-width", m.linewidth * 4)
                                    .attr("opacity", 0)
                                    .attr("class", "arrowline")
                                    .attr("ref", "arrow"+x1+"_"+x2+"_"+y1+"_"+y2)
                                    .on("click", make_editable)
                                    .on("dblclick", function(){
                                        var obj = d3.select(this);
                                        delete x.markers[d3.select(this.parentNode).attr("ref")];
                                        d3.select(this.parentNode).remove();
                                        self.svgDefault = false;
                                     });
                arrowline.start = arrowline.g.append("circle")
                                            .attr("r", 4)
                                            .attr("cx", xscale(x1)+margin.left)
                                            .attr("cy", y1)
                                            .attr("class", "arrowlineStart")
                                            .attr("opacity", 0)
                                            .call(d3.drag()
                                            .on("drag", function(d){
                                                var obj = d3.select(this);
                                                obj.style("cursor", "move");
                                                var coords = d3.mouse(svg.node());
                                                var posx = Math.round(xscale.invert(coords[0] - margin.left));
                                                var posy = Math.round(coords[1]);
                                                var parent = d3.select(this.parentNode);
                                                var ref = parent.attr("ref");
                                                var m = x.markers[ref];
                                                if(posx!=m.ref[0] || posy!=m.ref[1]){
                                                    m.ref = [posx, posy, m.ref[2], m.ref[3]];
                                                    x.markers["arrow"+posx+"_"+m.ref[2]+"_"+posy+"_"+m.ref[3]] = m;
                                                    delete(x.markers[ref]);
                                                    parent.attr("ref", "arrow"+posx+"_"+m.ref[2]+"_"+posy+"_"+m.ref[3]);
                                                    parent.selectAll(".arrowline")
                                                          .attr("x1", xscale(posx)+margin.left)
                                                          .attr("y1", posy)
                                                          .attr("ref", "arrow"+posx+"_"+m.ref[2]+"_"+posy+"_"+m.ref[3]);
                                                    obj.attr("cx", xscale(posx)+margin.left)
                                                       .attr("cy", posy);
                                                }
                                            })
                                            .on("end", function(d){
                                                 d3.select(this).style("cursor", "default");
                                            }));
                arrowline.end = arrowline.g.append("circle")
                                            .attr("r", 4)
                                            .attr("cx", xscale(x2)+margin.left)
                                            .attr("cy", y2)
                                            .attr("class", "arrowlineEnd")
                                            .attr("opacity", 0)
                                            .call(d3.drag()
                                            .on("drag", function(d){
                                                var obj = d3.select(this);
                                                obj.style("cursor", "move");
                                                var coords = d3.mouse(svg.node());
                                                var posx = Math.round(xscale.invert(coords[0] - margin.left));
                                                var posy = Math.round(coords[1]);
                                                var parent = d3.select(this.parentNode);
                                                var ref = parent.attr("ref");
                                                var m = x.markers[ref];
                                                if(posx!=m.ref[2] || posy!=m.ref[3]){
                                                    m.ref = [m.ref[0], m.ref[1], posx, posy];
                                                    x.markers["arrow"+m.ref[0]+"_"+posx+"_"+m.ref[1]+"_"+posy] = m;
                                                    delete(x.markers[ref]);
                                                    parent.attr("ref", "arrow"+m.ref[0]+"_"+posx+"_"+m.ref[1]+"_"+posy);
                                                    parent.selectAll(".arrowline")
                                                          .attr("x2", xscale(posx)+margin.left)
                                                          .attr("y2", posy)
                                                          .attr("ref", "arrow"+posx+"_"+m.ref[2]+"_"+posy+"_"+m.ref[3]);
                                                    obj.attr("cx", xscale(posx)+margin.left)
                                                       .attr("cy", posy);
                                                }
                                            })
                                            .on("end", function(d){
                                                 d3.select(this).style("cursor", "default");
                                            }));
                return(arrowline);
            };
            self.redraw = function(){
                self.g.remove();
                self.g = svg.insert("g");
                var rect = self.g.insert("rect").attr("width", svg.attr("width"))
                                     .attr("height", margin.top)
                                     .attr("fill", "white")
                                     .attr("stroke", "none")
                                     .attr("opacity", 0)
                                     .call(d3.drag()
                                             .on("start", self.dragstarted)
                                             .on("drag", self.dragged)
                                             .on("end", self.dragended));
                for(var k in x.markers){
                    var m = x.markers[k];
                    if(typeof(m)==="undefined") continue;
                    var l;
                    switch(m.markertype){
                        case 0:
                            l= self.g.append("line")
                                       .attr("y1", margin.top)
                                       .attr("y2", +svg.attr("height") - margin.bottom)
                                       .attr("x1", xscale(m.ref) + margin.left)
                                       .attr("x2", xscale(m.ref) + margin.left)
                                       .attr("stroke", m.color)
                                       .attr("stroke-width", "1px")
                                       .attr("fill", "none")
                                       .style("cursor", "move")
                                       .attr("ref", "line"+m.ref)
                                       .attr("class", "Marker");
                            if(m.linetype==="dashed"){
                                l.style("stroke-dasharray", ("3, 3"));
                            }else{
                                //solid
                            }
                            l.call(d3.drag().on("drag", self.draggedLine)
                                            .on("end", self.dragendLine));
                            break;
                        case 1:
                            l= self.g.append("rect")
                                   .attr("id", "rect"+m.ref)
                                   .attr("stroke", "none")
                                   .attr("fill", m.color)
                                   .attr("x", xscale(m.ref) + margin.left)
                                   .attr("width", xscale(m.ref+m.linewidth) - xscale(m.ref))
                                   .attr("y", margin.top)
                                   .attr("height", +svg.attr("height") - margin.bottom - margin.top)
                                   .attr("ref", "rect"+m.ref)
                                   .style("opacity", m.opacity)
                                   .attr("class", "Marker");
                            self.g.append("line")
                                   .attr("id", "rectL"+m.ref)
                                   .attr("y1", margin.top)
                                   .attr("y2", +svg.attr("height") - margin.bottom)
                                   .attr("x1", xscale(m.ref) + margin.left)
                                   .attr("x2", xscale(m.ref) + margin.left)
                                   .attr("stroke", "white")
                                   .attr("stroke-width", "2px")
                                   .style("opacity", 0)
                                   .style("cursor", "ew-resize")
                                   .attr("ref", "rect"+m.ref)
                                   .call(d3.drag().on("drag", self.resizeRectL));
                            self.g.append("line")
                                   .attr("id", "rectR"+m.ref)
                                   .attr("y1", margin.top)
                                   .attr("y2", +svg.attr("height") - margin.bottom)
                                   .attr("x1", xscale(m.ref+m.linewidth) + margin.left)
                                   .attr("x2", xscale(m.ref+m.linewidth) + margin.left)
                                   .attr("stroke", "white")
                                   .attr("stroke-width", "2px")
                                   .style("opacity", 0)
                                   .style("cursor", "ew-resize")
                                   .attr("ref", "rect"+m.ref)
                                   .call(d3.drag().on("drag", self.resizeRectR));
                            break;
                        case 2:
                            l= self.g.append("text")
                                   .attr("id", "text"+m.ref[0] + "_" + m.ref[1])
                                   .attr("ref", "text"+m.ref[0] + "_" + m.ref[1])
                                   .attr("fill", m.color)
                                   .attr("y", m.ref[1])
                                   .attr("x", xscale(m.ref[0]) + margin.left)
                                   .attr("text-anchor", "start")
                                   .style("font-size", m.fontsize+"px")
                                   .text(m.txt)
                                   .attr("class", "Mlabel")
                                   .call(makeLabelDraggable);
                            break;
                        case 3:
                            l= new Arrow(self.g, m);
                            break;
                    }
                    if(m.markertype!=3){
                        l.on("click", make_editable)
                         .on("dblclick", function(){
                            var obj = d3.select(this);
                            obj.remove();
                            delete x.markers[obj.attr("ref")];
                            self.redraw();
                            self.svgDefault = false;
                         });
                     }
                }
            }
            self.redraw();
            return(self);
        }
        
        var draw = function(vspace=10){
            if(typeof(g)!="undefined") g.remove();
            if(typeof(resizeBtn)!="undefined") resizeBtn.remove();
            if(typeof(marginBox)!="undefined") marginBox.remove();
            if(typeof(markerGroup)!="undefined") markerGroup.remove();
            xscale = d3.scaleLinear().domain([x.start, x.end]).rangeRound([0, widthF()]);
            g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
            g.append("rect")
             .attr("x", -margin.left)
             .attr("y", margin.top)
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
                //adjustable track height
                var trackBottom = track.append("line")
                                       .attr("kvalue", k)
                                       .attr("id", "trackbottom"+k)
                                       .attr("x1", 0)
                                       .attr("y1", thisHeight)
                                       .attr("x2", widthF())
                                       .attr("y2", thisHeight)
                                       .attr("stroke", "white")
                                       .attr("stroke-width", "3px")
                                       .attr("opacity", 0)
                                       .style("cursor", "ns-resize")
                                       .call(d3.drag()
                                       .on("start", function(d){
                                           coor = d3.mouse(svg.node());
                                       })
                                       .on("drag", function(d){
                                           var obj = d3.select(this);
                                           var dy = d3.mouse(svg.node())[1] - coor[1];
                                           coor = d3.mouse(svg.node());
                                           var k = +obj.attr("kvalue");
                                           var thisH = +obj.attr("y1") + dy;
                                           obj.attr("y1", thisH).attr("y2", thisH);
                                           var totalH = heightF() + dy;
                                           var ratio = heightF()/totalH;
                                           svg.attr("height", +svg.attr("height") + dy);
                                           for(var i=0; i<trackNames().length; i++){
                                               if(i==k){
                                                   x.height[trackNames()[i]] = (thisH + 2*vspace)/totalH;
                                               }else{
                                                   x.height[trackNames()[i]] = x.height[trackNames()[i]]*ratio;
                                               }
                                           }
                                           draw();
                                       }));
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
            d3.selectAll(".tick>text").on("click", make_editable);
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
            markerGroup = new Marker();
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
