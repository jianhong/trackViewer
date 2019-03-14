HTMLWidgets.widget({

  name: 'browseTracks',

  type: 'output',

  factory: function(el, width, height) {
  	// canvas
    var svg = d3.select(el).append("svg");
    // export menu
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
    var menu3 = d3.select(el).append("input")
                 .attr("type", "button")
                 .attr("name", "undo")
                 .attr("id", "undo")
                 .attr("value", "Undo");
    return {
      renderValue: function(x) {
        //export functions
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
            removeHighlight();
            
			if(typeof(cp)!="undefined"){
				cp.remove();
				cp = undefined;
			}
			
            svgAsDataUri(svg.node(), {width:svg.attr("width"), height:svg.attr("height"), scale: 4}, function(uri){
                var a = document.createElement('a');
                a.href = uri;
                a.download = 'trackViewer.svg';
                fireEvent(a, 'click');
            });
            plotregion.renew();
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
            removeHighlight();
            
			if(typeof(cp)!="undefined"){
				cp.remove();
				cp = undefined;
			}
			
            var background = svg.insert("rect", ":first-child")
                              .attr("width", svg.attr("width"))
                              .attr("height", svg.attr("height"))
                              .attr("fill", "white")
                              .attr("strock", "none");
            saveSvgAsPng(svg.node(), "trackViewer.png", {scale: 4});
            background.remove();
            plotregion.renew();
            ruler = new Ruler();
          });
        
        //set canvas size
        svg.attr("width", width)
           .attr("height", height);
        // set canvas margin
        var margin = {top: 20, right: 50, bottom: 30, left: 100};
        
        // default values
        var defaultFontSize = 16;
        var defaultColor = "#000000";
        
        //data x
        x.opacity = {};
        x.fontsize = {};//track y label font, track tick font, named as trackName__tick
        x.rotate = {};//track y label angle, track tick label angle, named as trackNames_tick
        x.color = {};//track y label color, track tick label color, named as trackNames_tick
        x.dataYlabelPos = {x:{},y:{}}; // track y label position
        x.lolliplotTrackLabel = {};
        
        //parameter to save all parameters
        var parameter = {};
        
        // get track name
        trackNames = function(){
        	if(typeof(x.name)=="string"){
        		x.name = [x.name];
        	}
            return(x.name);
        }; 
        safeNames = function(){
        	var n=clone(trackNames());
        	for(var i=0; i<n.length; i++){
        		n[i] = n[i].replace(/\.| /g, "_");
        	}
        	return(n);
        }
        
        for(var k=0; k<trackNames().length; k++){
            x.opacity[trackNames()[k]] = 1;
            x.fontsize[trackNames()[k]] = defaultFontSize;
        }
        
        var xHeight=function(){
            var xH=[];
            xH[0] = 0;
            for(var k=0; k<x.name.length; k++){
                xH[k+1] = x.height[trackNames()[k]]*heightF() + xH[k];
            }
            return(xH);
        };//get each track height position
                
        // select item
        var highlightItem;
        function removeHighlight(){
        	if(typeof(highlightItem)!="undefined"){
				highlightItem.rmHighlight();
			}
        }
        svg.on('click', removeHighlight);
        
        
         var changeTrackName = function(k, txt){
            if(x.name.indexOf(txt) != -1){
                console.log("used name " + txt);
                return;
            }
            var xkeys = d3.keys(x);
            var known = ["name", "chromosome", "start", "end", "strand", "markers"];
            xkeys = xkeys.filter(function(el){
                return known.indexOf(el)<0;
            });
            
            for(var i=0; i<xkeys.length; i++){
            	if(typeof(x[xkeys[i]][x.name[k]])!="undefined"){
            		x[xkeys[i]][txt] = x[xkeys[i]][x.name[k]];
					delete(x[xkeys[i]][x.name[k]]);
            	}
                
            }
            x.name[k] = txt;
        };
        
        function changeTranslate(ele, dx, dy){
			var s = ele.attr("transform");
			var t = s.substring(s.indexOf("(")+1, s.indexOf(")"))
					 .split(",");
			ele.attr("transform", function(){
				return("translate("+(Number(t[0])+dx)+","+(Number(t[1])+dy)+")");
			 });
        }
        //text
        var textDefaultOptions = function(){
        	return({
        	id: "textID", 
        	text: "label", 
        	x: 100, 
        	y: 100, 
        	vp: svg, 
        	dx: 0,
        	dy: 0,
        	angle: 0,
        	coor: 0,
        	trackKey: 0,
        	cls: "labelText",
        	fontsize: defaultFontSize,
        	color: "#000",
        	anchor: "middle",
        	datatrack: 0,
        	poskey:0,
        	colorPickerId:9,
        	ref: "text",
        	dblclick:function(){}
        	});
        
        };
        function Label(option=textDefaultOptions()){
            var self = this;
        	self.text = option.text;
        	self.x = option.x;
        	self.y = option.y;
        	self.dx = option.dx;
        	self.dy = option.dy;
        	self.oriX = option.x;
        	self.oriY = option.y;
        	self.anchor = option.anchor;
        	self.fontsize = option.fontsize;
        	self.color = option.color;
        	self.colorPickerId = option.colorPickerId;
        	self.vp = option.vp;
        	self.cls = option.cls;
        	self.k = option.trackKey;
        	self.coor = option.coor;
        	self.id = option.id;
        	self.angle = option.angle;
        	self.datatrack = option.datatrack;
        	self.poskey = option.poskey;
        	self.ref = option.ref;
        	self.dblclick = option.dblclick;
        	self.g = self.vp.append('g')
        	           .attr("id", self.id)
        	           .attr("transform", "translate("+self.x+","+self.y+")");
        	self.bgg = self.g.append('g').attr("transform", "rotate("+self.angle+")");
        	self.bg = self.bgg.append('g');
        	self.resize_bg = self.bgg.append('g');//resize point
        	self.rotate_bg = self.bgg.append('g');//rotate point
        	         
        	//move tracks order
        	var oriX,oriY;
			var tmp_id = []; // key: spring, value: parent
			self.dragstarted = function () {
			  d3.select(this).style("cursor", "move").raise().classed("active", true);
			  oriX = d3.event.x;
			  oriY = d3.event.y;
			  tmp_id.push({});
			  var tmpstatus = {
			  	ddx: parameter.xaxOpt.ddx,
			  	ddy: parameter.xaxOpt.ddy,
			  	ylabx0: x.dataYlabelPos.x[safeNames()[self.k]+"_0_"+self.ref+"_tick"],
			  	ylabx1: x.dataYlabelPos.x[safeNames()[self.k]+"_1_"+self.ref+"_tick"],
			  	ylaby0: x.dataYlabelPos.y[safeNames()[self.k]+"_0_"+self.ref+"_tick"],
			  	ylaby1: x.dataYlabelPos.y[safeNames()[self.k]+"_1_"+self.ref+"_tick"],
			  	yx: x.dataYlabelPos.x[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref],
			  	yy: x.dataYlabelPos.y[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref],
			  	label: clone(x.markers[self.id]),
			  	id: self.id,
			  	font: self.fontsize,
			  	lolliplotTrackLabel: x.lolliplotTrackLabel[safeNames()[self.k]+"_"+self.datatrack+"_"+self.poskey]
			  };
			  addNewHistory({
			  	undo:function(){
			  		if(self.cls=="xaxis_tick"){//xaxis tick
			  			parameter.xaxOpt.ddx = tmpstatus.ddx;
			  			parameter.xaxOpt.ddy = tmpstatus.ddy;
			  			parameter.xaxOpt.fontsize=tmpstatus.fontsize;
			  		}else{
			  			if(/yaxis_tick/.exec(self.cls)){//yaxis tick
							x.dataYlabelPos.x[safeNames()[self.k]+"_0_"+self.ref+"_tick"] = tmpstatus.ylabx0;
							x.dataYlabelPos.x[safeNames()[self.k]+"_1_"+self.ref+"_tick"] = tmpstatus.ylabx1;
							x.dataYlabelPos.y[safeNames()[self.k]+"_0_"+self.ref+"_tick"] = tmpstatus.ylaby0;
							x.dataYlabelPos.y[safeNames()[self.k]+"_1_"+self.ref+"_tick"] = tmpstatus.ylaby1;
							x.fontsize[safeNames()[self.k]+"__tick"]=tmpstatus.fontsize;
			  			}else{
			  				if(/dataYlabel_/.exec(self.cls)){
			  					x.dataYlabelPos.x[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref] = tmpstatus.yx;
			  					x.dataYlabelPos.y[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref] = tmpstatus.yy;
			  					x.fontsize[trackNames()[self.k]]=tmpstatus.fontsize;
			  				}else{
			  					if(/trackLayerLabel_/.exec(self.cls)){
			  						x.dataYlabelPos.x[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref] = tmpstatus.yx;
									x.dataYlabelPos.y[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref] = tmpstatus.yy;
									x.fontsize[safeNames()[self.k]+"_"+self.ref]=tmpstatus.fontsize;
			  					}else{
			  						if(self.cls=="Mlabel"){
			  							if(tmp_id.length>0){
			  							var keys = Object.keys(tmp_id[tmp_id.length-1]);
			  							for(var j=keys.length-1; j>=0; j--){
			  								if(tmp_id[tmp_id.length-1][keys[j]]==tmpstatus.id){
			  									if(typeof(x.markers[keys[j]])!="undefined"){
			  							 			delete(x.markers[keys[j]]);
												}
												break;
			  								}
			  							}
			  								tmp_id.pop();
			  							}
			  							x.markers[tmpstatus.id] = tmpstatus.label;
			  							x.markers[tmpstatus.id].fontsize=tmpstatus.font;
			  							
			  						}else{
			  							if(/nodelabel/.exec(self.cls)){
			  								x.lolliplotTrackLabel[safeNames()[self.k]+"_"+self.datatrack+"_"+self.poskey] = tmpstatus.lolliplotTrackLabel;
			  								x.fontsize["lolliplotTrackLabel_"+safeNames()[self.k]+"_"+self.datatrack+"_"+self.poskey]=tmpstatus.fontsize;
			  							}
			  						}
			  					}
			  				}
			  			}
			  		}
			  		plotregion.renew();
			  	},
			  	redo:function(){}
			  });
			};
			self.dragended = function () {
			  d3.select(this).style("cursor", "default").classed("active", false);
			};
			self.dragged = function() {
			  var dx = d3.event.x - oriX;
			  var dy = d3.event.y - oriY;
			  if(self.cls=="xaxis_tick"){//xaxis tick
					d3.selectAll(".xaxis_tick").each(function(){
						changeTranslate(d3.select(this.parentNode), dx, dy);
					});
					parameter.xaxOpt.ddx += dx;
					parameter.xaxOpt.ddy += dy;
				}else{                
					if(/yaxis_tick/.exec(self.cls)){//yaxis tick
						d3.selectAll("."+self.cls).each(function(){
							changeTranslate(d3.select(this.parentNode), dx, dy);
						});
						if(typeof(x.dataYlabelPos.x[safeNames()[self.k]+"_0_"+self.ref+"_tick"])=="undefined"){
							x.dataYlabelPos.x[safeNames()[self.k]+"_0_"+self.ref+"_tick"]=dx;
						}else{
							x.dataYlabelPos.x[safeNames()[self.k]+"_0_"+self.ref+"_tick"]+=dx;
						}
						if(typeof(x.dataYlabelPos.y[safeNames()[self.k]+"_0_"+self.ref+"_tick"])=="undefined"){
							x.dataYlabelPos.y[safeNames()[self.k]+"_0_"+self.ref+"_tick"]=dy;
						}else{
							x.dataYlabelPos.y[safeNames()[self.k]+"_0_"+self.ref+"_tick"]+=dy;
						}
						if(typeof(x.dataYlabelPos.x[safeNames()[self.k]+"_1_"+self.ref+"_tick"])=="undefined"){
							x.dataYlabelPos.x[safeNames()[self.k]+"_1_"+self.ref+"_tick"]=dx;
						}else{
							x.dataYlabelPos.x[safeNames()[self.k]+"_1_"+self.ref+"_tick"]+=dx;
						}
						if(typeof(x.dataYlabelPos.y[safeNames()[self.k]+"_1_"+self.ref+"_tick"])=="undefined"){
							x.dataYlabelPos.y[safeNames()[self.k]+"_1_"+self.ref+"_tick"]=dy;
						}else{
							x.dataYlabelPos.y[safeNames()[self.k]+"_1_"+self.ref+"_tick"]+=dy;
						}
					}else{
						  self.x=self.x+dx;
						  self.y=self.y+dy;
						if(/dataYlabel_/.exec(self.cls)){
							x.dataYlabelPos.x[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref]=xscale().invert(self.x);
							x.dataYlabelPos.y[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref]=self.y;
							changeTranslate(d3.select(this.parentNode), dx, dy);
						}else{
							if(/trackLayerLabel_/.exec(self.cls)){
								x.dataYlabelPos.x[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref]=xscale().invert(self.x);
								x.dataYlabelPos.y[safeNames()[self.k]+"_"+self.datatrack+"_"+self.ref]=self.y;
								changeTranslate(self.g, dx, dy);
							}else{
								if(self.cls=="Mlabel"){
									var coordinates = d3.mouse(svg.node());
									self.y = coordinates[1];
									self.x = coordinates[0];
									  var m = x.markers[self.id];
									  m.ref = [Math.round(xscale().invert(self.x-margin.left)), self.y];
									  delete(x.markers[self.id]);
									  var old_id = self.id;
									  self.id="text"+m.ref[0] + "_" + m.ref[1];
									  self.ref="text"+m.ref[0] + "_" + m.ref[1];
									  self.body.attr("ref", self.ref);
									  x.markers[self.id] = m;
									  if(typeof(tmp_id[tmp_id.length-1][old_id])!="undefined"){
									  	tmp_id[tmp_id.length-1][self.id] = tmp_id[tmp_id.length-1][old_id];
									  	if(old_id!=self.id) delete(tmp_id[tmp_id.length-1][old_id]);
									  }else{
									  	tmp_id[tmp_id.length-1][self.id] = old_id;
									  }
									  self.g.attr("transform", "translate("+self.x+","+self.y+")");
								}else{
									if(/nodelabel/.exec(self.cls)){
										if(typeof(x.lolliplotTrackLabel[safeNames()[self.k]+"_"+self.datatrack+"_"+self.poskey])!="undefined"){
											x.lolliplotTrackLabel[safeNames()[self.k]+"_"+self.datatrack+"_"+self.poskey].x += dx;
											x.lolliplotTrackLabel[safeNames()[self.k]+"_"+self.datatrack+"_"+self.poskey].y += dy;
										}else{
											x.lolliplotTrackLabel[safeNames()[self.k]+"_"+self.datatrack+"_"+self.poskey] = {x:dx, y:dy};
										}
										changeTranslate(d3.select(this.parentNode), dx, dy);
									}else{
										changeTranslate(d3.select(this.parentNode), dx, dy);
									}
								}
							}
						}
					}
				}
				
			};
			self.resize_draged = function(){
                self.fontsize = self.fontsize - d3.event.dy;
                if(self.cls=="xaxis_tick"){//xaxis tick
                	d3.selectAll(".xaxis_tick").style("font-size", self.fontsize+"px");
                	parameter.xaxOpt.fontsize=self.fontsize;
                }else{                
                	if(/yaxis_tick/.exec(self.cls)){//yaxis tick
                		d3.selectAll("."+self.cls).style("font-size", self.fontsize+"px");
                		x.fontsize[safeNames()[self.k]+"__tick"]=self.fontsize;
                	}else{
                		if(/dataYlabel_/.exec(self.cls)){
                			x.fontsize[trackNames()[self.k]]=self.fontsize;
                		}
                		if(/trackLayerLabel/.exec(self.cls)){
                			x.fontsize[safeNames()[self.k]+"_"+self.ref]=self.fontsize;
                		}
                		if(/nodelabel/.exec(self.cls)){
                			x.fontsize["lolliplotTrackLabel_"+safeNames()[self.k]+"_"+self.datatrack+"_"+self.poskey]=self.fontsize;
                		}
                		if(/Mlabel/.exec(self.cls)){
                			x.markers[self.id].fontsize=self.fontsize;
                		}
                		self.body.style("font-size", self.fontsize+"px");
                	}
                }
                var SVGRect = self.body.node().getBBox();
                self.bg.attr("x", SVGRect.x)
                        .attr("y", SVGRect.y)
                        .attr("width", SVGRect.width)
                        .attr("height", SVGRect.height);
                self.resize_bg.attr("transform", "translate("+(SVGRect.x + SVGRect.width/2)+","+SVGRect.y+")");
                self.rotate_bg.attr("transform", "translate("+(SVGRect.x + SVGRect.width)+","+SVGRect.y+")");
            };
			self.onClick = function (){
				removeHighlight();
				self.highlight();
				self.editable();
				if(!/axis_tick/.exec(self.cls)){ColorPicker(this, self.colorPickerId);}
			};
						
			self.rmHighlight = function(){
				self.bg.remove();
			    self.resize_bg.remove();
			    self.rotate_bg.remove();
			};
			self.highlight = function(){
			  if(typeof(self.g)!="undefined"){
			    self.rmHighlight();
			  	var SVGRect = self.body.node().getBBox();
				self.bg = self.bgg.insert("rect", ":first-child")
                        .attr("x", SVGRect.x)
                        .attr("y", SVGRect.y)
                        .attr("width", SVGRect.width)
                        .attr("height", SVGRect.height)
                        .attr("fill", "#B3D8FD");
                self.resizeable();
                self.rotatable();
                highlightItem=self;
              }
            };
            
			self.resizeable = function(){
			  if(typeof(self.g)!="undefined"){
			  	var SVGRect = self.body.node().getBBox();
				var handlerRadius = 2.5;
                self.resize_bg = self.bgg.append("circle")
                          .style("fill", "white")
                          .style("stroke", "blue")
                          .attr("transform", "translate("+(SVGRect.x + SVGRect.width/2)+","+SVGRect.y+")")
                          .attr("r", handlerRadius)
                          .call(d3.drag()
                                      .on("start", self.dragstarted)
                                      .on("drag", self.resize_draged)
                                      .on("end", self.dragended));
              }
			};
			function rotateText(ang){
				var tmpstatus = self.angle;
				addNewHistory({
					undo:function(){
						self.angle = tmpstatus;
						if(self.cls=="xaxis_tick"){//xaxis tick
							parameter.xaxOpt.angle=tmpstatus;
						}else{
							if(/yaxis_tick/.exec(self.cls)){//yaxis tick
								x.rotate[safeNames()[self.k]+"__tick"] = tmpstatus;
							}else{
								if(/dataYlabel_/.exec(self.cls)){
									x.rotate[trackNames()[self.k]]=tmpstatus;
								}
								if(/trackLayerLabel/.exec(self.cls)){
									x.rotate[safeNames()[self.k]+"_"+self.ref]=tmpstatus;
								}
								if(/nodelabel/.exec(self.cls)){
									x.tracklist[trackNames()[self.k]][self.datatrack]["label.parameter.rot"][self.poskey]=-tmpstatus;
								}
								if(/Mlabel/.exec(self.cls)){
									x.markers[self.id].angle=tmpstatus;
								}
							}
						}
						plotregion.renew();
					},
					redo:function(){}
				});
				self.angle = self.angle + ang;
				if(self.cls=="xaxis_tick"){//xaxis tick
                	d3.selectAll(".xaxis_tick")
                	       .attr("transform", "rotate("+self.angle+")");
                	parameter.xaxOpt.angle=self.angle;
                }else{                
                	if(/yaxis_tick/.exec(self.cls)){//yaxis tick
                		d3.selectAll("."+self.cls).attr("transform", "rotate("+self.angle+")");
                		x.rotate[safeNames()[self.k]+"__tick"] = self.angle;
                	}else{
						if(/dataYlabel_/.exec(self.cls)){
                			x.rotate[trackNames()[self.k]]=self.angle;
                		}
                		if(/trackLayerLabel/.exec(self.cls)){
                			x.rotate[safeNames()[self.k]+"_"+self.ref]=self.angle;
                		}
                		if(/nodelabel/.exec(self.cls)){
                			if(typeof(x.tracklist[trackNames()[self.k]][self.datatrack]["label.parameter.rot"])=="undefined"){
                				x.tracklist[trackNames()[self.k]][self.datatrack]["label.parameter.rot"] = [];
                				for(var i=0; i<x.tracklist[trackNames()[self.k]][self.datatrack].start.length; i++){
                					x.tracklist[trackNames()[self.k]][self.datatrack]["label.parameter.rot"].push(90);
                				}
                			}
                			x.tracklist[trackNames()[self.k]][self.datatrack]["label.parameter.rot"][self.poskey]=-self.angle;
                		}
                		if(/Mlabel/.exec(self.cls)){
                			x.markers[self.id].angle=self.angle;
                		}
                	}
                }
				self.body.attr("transform", "rotate("+self.angle+")");
				self.bgg.attr("transform", "rotate("+self.angle+")");
			}
			self.rotatable = function(){
			  if(typeof(self.g)!="undefined"){
			  	var SVGRect = self.body.node().getBBox();
                self.rotate_bg = self.bgg.append('g')
                                       .attr("transform", "translate("+(SVGRect.x + SVGRect.width)+","+SVGRect.y+")");
                self.rotate_bg.append("path")
						.attr('d', function(d){
							return 'M 0 0 l 10 -10 l -10 0 z';
						})
						.style('fill', "black")
						.style('stroke', 'black')
						.on("click", function(){
							rotateText(-45);
						});
				self.rotate_bg.append("path")
						.attr('d', function(d){
							return 'M 0 0 l 10 -10 l 0 10 z';
						})
						.style('fill', "white")
						.style('stroke', 'black')
						.on("click", function(){
							rotateText(45);
						});
              }
			};
			self.editable = function(){
				// inject a HTML form to edit the content here...
				var xy = self.body.node().getBBox();
				var p_xy = self.g.node().getBBox();

				xy.x -= p_xy.x/2;
				xy.y -= p_xy.y/2;

				var frm = self.g.append("foreignObject");
								
				var inp = frm
					.attr("x", xy.x)
					.attr("y", xy.y)
					.attr("width", 300)
					.attr("height", 25)
					.append("xhtml:form")
						.append("input")
							.attr("value", function() {
								this.focus();
								return self.text;
							})
							.attr("style", "width: 200px;")
							// make the form go away when you jump out (form looses focus) or hit ENTER:
							.on("blur", function() {
								var tmpstatus = clone(self);
								if(self.text != inp.node().value){
									addNewHistory({
										undo:function(){
											if(/yaxis_/.exec(tmpstatus.cls)){
												x.tracklist[trackNames()[tmpstatus.k]].ylim[tmpstatus.datatrack]=Number(tmpstatus.text);
											}
											if(/dataYlabel_/.exec(tmpstatus.cls)){
												changeTrackName(tmpstatus.k, tmpstatus.text)
											}
											if(/Mlabel/.exec(tmpstatus.cls)){
												  x.markers[tmpstatus.id].ref = [xscale().invert(tmpstatus.x-margin.left), tmpstatus.y];
												  x.markers[tmpstatus.id].text = tmpstatus.text;
											}
											if(/legend/.exec(tmpstatus.cls)){
												x.markers[tmpstatus.id].text = tmpstatus.text;
											}
											if(/trackLayerLabel/.exec(tmpstatus.cls)){
												parameter.trackLayerDataTxt[trackNames()[tmpstatus.k]+"_"+tmpstatus.ref] = tmpstatus.text;
											}
											plotregion.renew();
										},
										redo:function(){}
									});
								}
								var txt = inp.node().value;
								var old = self.text;
								self.text = txt;
								self.body.text(self.text);
								self.highlight();
								frm.remove();
								if(/yaxis_/.exec(self.cls)){
									x.tracklist[trackNames()[self.k]].ylim[self.datatrack]=Number(self.text);
									if(Math.abs(old-self.text)>2) plotregion.renew();
								}
								if(/dataYlabel_/.exec(self.cls)){
									changeTrackName(self.k, self.text)
								}
								if(/Mlabel/.exec(self.cls)){
									  x.markers[self.id].text = txt;  
								}
								if(/legend/.exec(self.cls)){
									x.markers[self.id].text = txt;
								}
								if(/trackLayerLabel/.exec(self.cls)){
									parameter.trackLayerDataTxt[safeNames()[self.k]+"_"+self.ref] = txt;
								}
							})
							.on("keypress", function() {
								// IE fix
								if (!d3.event)
									d3.event = window.event;

								var e = d3.event;
								if (e.keyCode == 13)
								{
									if (typeof(e.cancelBubble) !== 'undefined') // IE
									  e.cancelBubble = true;
									if (e.stopPropagation)
									  e.stopPropagation();
									e.preventDefault();
									
									
									var tmpstatus = clone(self);
									if(self.text != inp.node().value){
										addNewHistory({
											undo:function(){
												if(/yaxis_/.exec(tmpstatus.cls)){
													x.tracklist[trackNames()[tmpstatus.k]].ylim[tmpstatus.datatrack]=Number(tmpstatus.text);
												}
												if(/dataYlabel_/.exec(tmpstatus.cls)){
													changeTrackName(tmpstatus.k, tmpstatus.text)
												}
												if(/Mlabel/.exec(tmpstatus.cls)){
													  x.markers[tmpstatus.id].ref = [xscale().invert(tmpstatus.x-margin.left), tmpstatus.y];
													  x.markers[tmpstatus.id].text = tmpstatus.text;
												}
												if(/legend/.exec(tmpstatus.cls)){
													x.markers[tmpstatus.id].text = tmpstatus.text;
												}
												if(/trackLayerLabel/.exec(tmpstatus.cls)){
													parameter.trackLayerDataTxt[trackNames()[tmpstatus.k]+"_"+tmpstatus.ref] = tmpstatus.text;
												}
												plotregion.renew();
											},
											redo:function(){}
										});
									}
									
									var old = self.text;
									self.text = inp.node().value;
									self.body.text(self.text);
									self.highlight();
									try{
										frm.remove();
									}catch(err){
									   //do nothing.
									}
									
									if(/yaxis_/.exec(self.cls)){
										x.tracklist[trackNames()[self.k]].ylim[self.datatrack]=Number(self.text);
										if(Math.abs(old-self.text)>2) plotregion.renew();
									}
									if(/dataYlabel_/.exec(self.cls)){
										changeTrackName(self.k, self.text)
									}
									if(/Mlabel/.exec(self.cls)){
										x.markers[self.id].text = self.text;
									}
									if(/legend/.exec(self.cls)){
										x.markers[self.id].text = self.text;
									}	
									if(/trackLayerLabel/.exec(self.cls)){
										parameter.trackLayerDataTxt[safeNames()[self.k]+"_"+self.ref] = self.text;
									}		
								}
							});
			};
			
        	self.body = self.g.append("text")
        		  .attr("x", 0)
        		  .attr("y", 0)
        		  .attr("dx", self.dx)
        		  .attr("dy", self.dy)
        		  .attr("text-anchor", self.anchor)
        		  .attr("class", self.cls)
        		  .attr("kvalue", self.k)
        		  .attr("datatrack", self.datatrack)
        		  .attr("poskey", self.poskey)
        		  .attr("fill", self.color)
        		  .attr("ref", self.ref)
        		  .text(self.text)
        		  .style("font-size", self.fontsize+"px")
        		  .attr("transform", "rotate("+self.angle+")")
        		  .on("click", self.onClick)
        		  .on("dblclick", self.dblclick)
        		  .call(d3.drag().on("start", self.dragstarted)
								 .on("drag", self.dragged)
								 .on("end", self.dragended));
			self.remove = function(){
				self.body.remove();
				self.g.remove();
			};
			
			return(self);
        }
        
        var legendLabel = function(option=textDefaultOptions()){
        	var self=this;
        	self.text = option.text;
        	self.x = option.x;
        	self.y = option.y;
        	self.dx = option.dx;
        	self.dy = option.dy;
        	self.anchor = option.anchor;
        	self.fontsize = option.fontsize;
        	self.color = option.color;
        	self.colorPickerId = option.colorPickerId;
        	self.cls = option.cls;
        	self.trackKey = option.trackKey || currentLayer;
        	self.coor = option.coor;
        	self.id = option.id;
        	self.angle = option.angle;
        	self.datatrack = option.datatrack;
        	self.poskey = option.poskey;
        	self.ref = option.ref;
        	self.fill = option.fill || "black";
        	self.border = option.border || "black";
        	self.r = option.r || 8;
        	self.vp = option.vp || svg;
        	
        	var M = clone(x.markers[self.id]);
        	var oriX,oriY,curX,curY;
			self.vpg = self.vp.append("g")
						 .attr("transform", "translate("+xscale()(self.x)+","+self.y+")")
						 .attr("id", self.id)
						 .call(d3.drag().on("drag", function(d){
								var obj = d3.select(this);
								obj.style("cursor", "move");
								var coords = d3.mouse(svg.node());
								var dx = coords[0] - curX;
								var dy = coords[1] - curY;
								curX = coords[0];
								curY = coords[1];
								changeTranslate(obj, dx, dy);
						 }).on("start", function(d){
						 	d3.select(this).style("cursor", "move").raise().classed("active", true);
						 	var coords = d3.mouse(svg.node());
						 	oriX = coords[0];
						 	oriY = coords[1];
							curX = coords[0];
							curY = coords[1];
						 	var tmpstatus={};
						 	if(typeof(x.markers[M.id])!="undefined"){
								tmpstatus={
									x: x.markers[M.id].x,
									y: x.markers[M.id].y
								};
						 	}
							addNewHistory({
								undo:function(){
										if(typeof(x.markers[M.id])!="undefined"){
											x.markers[M.id].x = tmpstatus.x;
											x.markers[M.id].y = tmpstatus.y;
											plotregion.renew();
										}
									},
								redo:function(){}
							});
						 }).on("end", function(d){
						 	d3.select(this).style("cursor", "default").classed("active", false);
						 	var m = x.markers[d3.select(this).attr("id")];
						 	m.x = Math.round(xscale().invert(curX));
						 	m.y = Math.round(curY);
						 	x.markers[d3.select(this).attr("id")] = m;
						 }))
						 .on("dblclick", self.remove);
	 
        	self.remove = function(){	
        		addNewHistory({
        			undo: function(){
        				x.markers[M.id] = M;
						plotregion.renew();
        			},
        			redo: self.remove
        		});
        		delete x.markers[self.ref];
        		self.vpg.remove();
        	};
        	
        	var marker = self;
        	marker.x =0;
        	marker.y =0;
        	marker.vp = self.vpg;
        	self.marker = new Label(marker);
        	self.circle = self.vpg.append("circle")
                				.attr("kvalue", self.trackKey)
                				.attr("cx", -self.fontsize/2-2)
                				.attr("cy", -self.fontsize/3.5)
                				.attr("r", self.r)
                				.attr("fill", self.fill)
                				.attr("stroke", self.border)
                				.attr("ref", self.ref)
                				.on("click", function(){
                					var obj = d3.select(this);
                					var picked = function(col){
										var tmpstatus = {k:self.id, v:x.markers[self.id].fill};
										addNewHistory({
											undo:function(){
												x.markers[tmpstatus.k].fill = tmpstatus.v;
												plotregion.renew();
											},
											redo:function(){}
										});
                						var m = x.markers[self.id];
                						m.fill = col;
                						x.markers[self.id] = m;
                						obj.attr("fill", col);
                					};
                					ColorPicker(this, picked);
                				})
						 .on("dblclick", self.remove);

			return(self);
		};
		
        //Markers
        x.markers = [];
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
            /*var makeLabelDraggable = d3.drag()
                       .on("drag", function(d){
                            var curr = d3.select(this);
                            if(typeof(curr.attr("ref"))!="undefined"){
                              var coords = d3.mouse(svg.node());
                              var posx = Math.round(xscale().invert(coords[0]-margin.left));
                              var posy = Math.round(coords[1]);
                              var m = x.markers[curr.attr("ref")];
                              if(posx!=m.ref[0] && posy!=m.ref[1]){
                                  m.ref = [posx, posy];
                                  x.markers["text"+posx + "_" + posy] = m;
                                  delete(x.markers[curr.attr("ref")]);
                                  curr.attr("id", "text"+m.ref[0] + "_" + m.ref[1])
                                     .attr("ref", "text"+m.ref[0] + "_" + m.ref[1])
                                     .attr("x", xscale()(m.ref[0]) + margin.left)
                                     .attr("y", m.ref[1])
                                     .style("cursor", "move");
                              }
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
                       });*/
                       
            self.addArrowLabel = function(x0, y0){
              // react on right-clicking
                if(self.svgDefault){
                    var coords = [x0, y0];//d3.mouse(svg.node());
                    var posx = Math.round(xscale().invert(coords[0]-margin.left+5));
                    var posy = Math.round(coords[1])-8; // can not keep position
                    var posx1 = Math.round(xscale().invert(coords[0]-margin.left+20));
                    var posy1 = posy -20;
                    var m = {
                        "ref" : [posx1, posy1],
                        "markertype" : 2,
                        "color" : "black",
                        "opacity" : 1,
                        "fontsize" : 12,
                        "text" : "label",
                        "angle": 0
                    };
                    var m1 = {
                        "ref" : [posx1, posy1, posx, posy],
                        "markertype" : 3,
                        "color" : "#000",
                        "opacity" : 1,
                        "linetype" : "solid",
                        "linewidth" : 1
                    };
                    
                	var tmpstatus={
                		text: "text"+posx1 + "_" + posy1,
                		arrow: "arrow"+m1.ref[0]+"_"+m1.ref[2]+"_"+m1.ref[1]+"_"+m1.ref[3]
                	};
                	addNewHistory({
                		undo:function(){
                			if(typeof(x.markers[tmpstatus.text])!="undefined") delete x.markers[tmpstatus.text];
                			if(typeof(x.markers[tmpstatus.arrow])!="undefined") delete x.markers[tmpstatus.arrow];
                			plotregion.renew();
                		},
                		redo:function(){}
                	});
                    x.markers["text"+posx1 + "_" + posy1] = m;
                    var opt = textDefaultOptions();
                    opt.id = "text"+m.ref[0] + "_" + m.ref[1];
                    opt.ref = "text"+m.ref[0] + "_" + m.ref[1];
                    opt.color = m.color;
                    opt.x = xscale()(m.ref[0]) + margin.left;
                    opt.y = m.ref[1];
                    opt.anchor = "start";
                    opt.fontsize = m.fontsize;
                    opt.angle = m.angle;
                    opt.text = m.text;
                    opt.vp = self.g;
                    opt.cls = "Mlabel";
                    opt.colorPickerId = 9;
                    var t = new Label(opt);
                    x.markers["arrow"+m1.ref[0]+"_"+m1.ref[2]+"_"+m1.ref[1]+"_"+m1.ref[3]]=m1;
                    var arrow = new Arrow(self.g, m1);
                }else{
                    self.svgDefault=true;
                }
            };
            self.addLegendLabel = function(x0, y0){
            	if(self.svgDefault){
                    var coords = [x0, y0];//d3.mouse(svg.node());
                    var m = textDefaultOptions();
                    m.text = "legend";
                    m.x = Math.round(xscale().invert(coords[0]));
                    m.y = Math.round(coords[1]);
                    m.anchor = "start";
                    m.color = "black";
                    m.cls = "legend_"+currentLayer+"_"+0;
					m.trackKey = k;
					m.id = "legend_"+currentLayer+"_"+0+"_"+m.color+"_"+m.text;
					m.angle = 0;
					m.datatrack = 0;
					m.ref = m.id;
					m.fill = m.color;
					m.border = "black";
					m.r = m.fontsize/2;
					m.vp = self.g;
					m.markertype = 4;
                    x.markers[m.id]=m;
                    var tmpstatus=m.id;
                	addNewHistory({
                		undo:function(){
                			if(typeof(x.markers[tmpstatus])!="undefined") delete x.markers[tmpstatus];
                			plotregion.renew();
                		},
                		redo:function(){}
                	});
                	
                    var llab = new legendLabel(m);
                }else{
                    self.svgDefault=true;
                }
            };
            
            self.g = svg.append("g");
            self.draggedLine = function(d){
                var coords = d3.mouse(svg.node())[0];
                d3.select(this).attr("x1", coords).attr("x2", coords);
                d3.select("line[ref='"+d3.select(this).attr("ref")+"']")
                	.attr("x1", coords).attr("x2", coords);
            }
            self.dragendLine = function(d){
                var coords = d3.mouse(svg.node());
                var posx = Math.round(xscale().invert(coords[0]-margin.left));
                var old = d3.select(this).attr("ref");
                var oldx = Number(old.replace("line", ""));
                if(Math.abs(xscale()(oldx)+ margin.left - coords[0])>1){
                    x.markers["line"+posx] = x.markers[old];
                    delete(x.markers[old]);
                    x.markers["line"+posx].ref = posx;
                    d3.select(this).attr("ref", "line"+posx);
                    self.redraw();
                    addNewHistory({
                		undo:function(){
                			x.markers[old] = x.markers["line"+posx];
                			x.markers["line"+posx].ref = Number(old.replace("line", ""));
                			delete(x.markers["line"+posx]);
                			self.redraw();
                		},
                		redo:function(){}
                	});
                }else{
                	var obj = d3.select(this);
					var picked = function(col){
						var tmpstatus = {k:obj.attr("ref"), v:x.markers[obj.attr("ref")].color};
						addNewHistory({
							undo:function(){
								x.markers[tmpstatus.k].color = tmpstatus.v;
								plotregion.renew();
							},
							redo:function(){}
						});
						x.markers[obj.attr("ref")].color = col;
						self.redraw();
					};
					ColorPicker(this, picked);
                }
            }
            self.dragstarted = function(d){
                            coor = d3.mouse(this);
                            var posx = Math.round(xscale().invert(coor[0] - margin.left));
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
                                           .attr("x", xscale()(m.ref) + margin.left)
                                           .attr("width", m.linewidth)
                                           .attr("y", margin.top)
                                           .attr("height", +svg.attr("height") - margin.bottom - margin.top)
                                           .attr("ref", "rect"+m.ref)
                                           .style("opacity", m.opacity)
                                           .attr("class", "Marker");
                };
              self.dragged = function(d){
                                    var coords = d3.mouse(this);
                                    var posx = Math.round(xscale().invert(coor[0] - margin.left));
                                    d3.select("#rect"+posx).attr("width", coords[0] - coor[0]);
                                    };
              self.dragended = function(d){
                                var coords = d3.mouse(this);
                                var posx = Math.round(xscale().invert(coor[0] - margin.left));
                                if(coords[0] - coor[0] > 3){
                                    x.markers["rect"+posx] = {
                                        "ref" : posx,
                                        "markertype" : 1,
                                        "color" : "gray",
                                        "opacity" : 0.3,
                                        "linetype" : "solid",
                                        "linewidth" : xscale().invert(coords[0] - margin.left) - posx
                                    };
                                    addNewHistory({
										undo:function(){
											delete(x.markers["rect"+posx]);
											self.redraw();
										},
										redo:function(){}
									});
                                }else{
                                    x.markers["line"+posx] = {
                                        "ref" : posx,
                                        "markertype" : 0,
                                        "color" : "gray",
                                        "opacity" : 1,
                                        "linetype" : "dashed",
                                        "linewidth" : 1
                                    };
                                    addNewHistory({
										undo:function(){
											delete(x.markers["line"+posx]);
											self.redraw();
										},
										redo:function(){}
									});
                                }
                                self.redraw();
                            };
            self.resizeRectL = function(d){
                var coords = d3.mouse(this);
                var posx = Math.round(xscale().invert(coords[0] - margin.left));
                var obj = d3.select(this);
                var ref = obj.attr("ref");
                if(posx != x.markers[ref].ref){
                    obj.attr("x1", xscale()(posx) + margin.left)
                       .attr("x2", xscale()(posx) + margin.left);
                    var newWid = d3.select("#rectR"+x.markers[ref].ref).attr("x1") - coords[0];
                    d3.select("#rect"+x.markers[ref].ref)
                      .attr("x", xscale()(posx) + margin.left)
                      .attr("width", newWid>0?newWid:0);
                }
            };
            self.resizeRectLEnd = function(d){
                var coords = d3.mouse(this);
                var posx = Math.round(xscale().invert(coords[0] - margin.left));
                var obj = d3.select(this);
                var ref = obj.attr("ref");
                if(posx != x.markers[ref].ref){
                	var tmpstatus = {k:ref,k2:posx,v:clone(x.markers[ref])};
                	addNewHistory({
                		undo:function(){
                			x.markers[tmpstatus.k] = clone(tmpstatus.v);
                			delete(x.markers["rect"+tmpstatus.k2]);
                			plotregion.renew();
                		},
                		redo:function(){}
                	});
                    x.markers[ref].linewidth = x.markers[ref].ref - posx + x.markers[ref].linewidth;
                    obj.attr("ref", "rect" + posx)
                       .attr("id", "rectL" + posx)
                       .attr("x1", xscale()(posx) + margin.left)
                       .attr("x2", xscale()(posx) + margin.left);
                    var newWid = xscale()(posx + x.markers[ref].linewidth) - xscale()(posx);
                    d3.select("#rect"+x.markers[ref].ref)
                      .attr("ref", "rect" + posx)
                      .attr("id", "rect" + posx)
                      .attr("x", xscale()(posx) + margin.left)
                      .attr("width", newWid>0?newWid:0);
                    d3.select("#rectR"+x.markers[ref].ref)
                      .attr("ref", "rect" + posx)
                      .attr("id", "rectR" + posx)
                      .attr("x1", xscale()(x.markers[ref].linewidth + posx)+margin.left)
                      .attr("x2", xscale()(x.markers[ref].linewidth + posx)+margin.left)
                    x.markers[ref].ref = posx;
                    x.markers["rect"+posx] = x.markers[ref];
                    delete(x.markers[ref]);
                    if(newWid<=0){
                      obj.remove();
                    }
                }
            };
            self.resizeRectR = function(d){
                var coords = d3.mouse(this);
                var posx = Math.round(xscale().invert(coords[0] - margin.left));
                var obj = d3.select(this);
                var ref = obj.attr("ref");
				obj.attr("x1", xscale()(posx) + margin.left)
				   .attr("x2", xscale()(posx) + margin.left);
				var newWid = xscale()(posx) - xscale()(x.markers[ref].ref);
				d3.select("#rect"+x.markers[ref].ref)
				  .attr("width", newWid>0?newWid:0);
            };
            self.resizeRectREnd = function(d){
                var coords = d3.mouse(this);
                var posx = Math.round(xscale().invert(coords[0] - margin.left));
                var obj = d3.select(this);
                var ref = obj.attr("ref");
                if(posx != x.markers[ref].ref + x.markers[ref].linewidth){
                	var tmpstatus = {k:ref,v:clone(x.markers[ref])};
                	addNewHistory({
                		undo:function(){
                			x.markers[tmpstatus.k] = clone(tmpstatus.v);
                			plotregion.renew();
                		},
                		redo:function(){}
                	});
                    obj.attr("x1", xscale()(posx) + margin.left)
                       .attr("x2", xscale()(posx) + margin.left);
                    var newWid = xscale()(posx) - xscale()(x.markers[ref].ref);
                    d3.select("#rect"+x.markers[ref].ref)
                      .attr("width", newWid>0?newWid:0);
                    x.markers[ref].linewidth = posx - x.markers[ref].ref;
                    if(newWid<=0){
                      obj.remove();
                    }
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
                                    .attr("x1", xscale()(x1)+margin.left)
                                    .attr("x2", xscale()(x2)+margin.left)
                                    .attr("y1", y1)
                                    .attr("y2", y2)
                                    .attr("stroke", m.color)
                                    .attr("stroke-width", m.linewidth)
                                    .attr("opacity", m.opacity)
                                    .attr("class", "arrowline")
                                    .attr("ref", "arrow"+x1+"_"+x2+"_"+y1+"_"+y2)
                                    .attr("marker-end", "url(#arrowhead"+x1+"_"+x2+"_"+y1+"_"+y2+")");
                arrowline.linecover = arrowline.g.append("line")
                                    .attr("x1", xscale()(x1)+margin.left)
                                    .attr("x2", xscale()(x2)+margin.left)
                                    .attr("y1", y1)
                                    .attr("y2", y2)
                                    .attr("stroke", m.color)
                                    .attr("stroke-width", m.linewidth * 10)
                                    .attr("opacity", 0)
                                    .attr("class", "arrowline")
                                    .attr("ref", "arrow"+x1+"_"+x2+"_"+y1+"_"+y2)
                                    .on("click", function(){
                                    	var obj = d3.select(this);
                                    	var picked = function(col){
											var tmpstatus = {k:obj.attr("ref"), v:x.markers[obj.attr("ref")].color};
											addNewHistory({
												undo:function(){
													x.markers[tmpstatus.k].color = tmpstatus.v;
													plotregion.renew();
												},
												redo:function(){}
											});
                                    		x.markers[obj.attr("ref")].color = col;
                                    		self.redraw();
                                    	};
                                    	ColorPicker(this, picked);
                                    })
                                    .on("dblclick", function(){
                                        var obj = d3.select(this);
                                        var par = d3.select(this.parentNode);
                                        var ref = par.attr("ref");
                                    	addNewHistory({
                                    		undo:function(){
                                    			var m1 = ref.replace("arrow", "").split("_");
												var m={
													"ref" : [m1[0], m1[2], m1[1], m1[3]],
													"markertype" : 3,
													"color" : "#000",
													"opacity" : 1,
													"linetype" : "solid",
													"linewidth" : 1
												};
                                    			x.markers[ref]=m;
                                    			self.redraw();
                                    		},
                                    		redo:function(){}
                                    	});
                                        delete x.markers[ref];
                                        par.remove();
                                        self.svgDefault = false;
                                     });
                var oldstart={x1:x1, y1:y1, x2:x2, y2:y2, col:m.color};
                arrowline.start = arrowline.g.append("circle")
                                            .attr("r", 10)
                                            .attr("cx", xscale()(x1)+margin.left)
                                            .attr("cy", y1)
                                            .attr("class", "arrowlineStart")
                                            .attr("opacity", 0)
                                            .call(d3.drag()
                                            .on("drag", function(d){
                                                var obj = d3.select(this);
                                                obj.style("cursor", "move");
                                                var coords = d3.mouse(svg.node());
                                                var posx = Math.round(xscale().invert(coords[0] - margin.left));
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
                                                          .attr("x1", xscale()(posx)+margin.left)
                                                          .attr("y1", posy)
                                                          .attr("ref", "arrow"+posx+"_"+m.ref[2]+"_"+posy+"_"+m.ref[3]);
                                                    obj.attr("cx", xscale()(posx)+margin.left)
                                                       .attr("cy", posy);
                                                }
                                            })
                                            .on("end", function(d){
                                                 d3.select(this).style("cursor", "default");
                                                 var ref=d3.select(this.parentNode).attr("ref");
                                                 addNewHistory({
													undo:function(){
														delete(x.markers[ref]);
														var m={
															"ref" : [oldstart.x1, oldstart.y1, oldstart.x2, oldstart.y2],
															"markertype" : 3,
															"color" : oldstart.col,
															"opacity" : 1,
															"linetype" : "solid",
															"linewidth" : 1
														};
														x.markers["arrow"+oldstart.x1+"_"+oldstart.x2+"_"+oldstart.y1+"_"+oldstart.y2]=m;
														self.redraw();
													},
													redo:function(){}
												});
                                            })
                                            .on("start", function(d){
                                            	d3.select(this).style("cursor", "move");
                                            	oldstart={x1:x1, y1:y1, x2:x2, y2:y2, col:m.color};
                                            }));
                arrowline.end = arrowline.g.append("circle")
                                            .attr("r", 10)
                                            .attr("cx", xscale()(x2)+margin.left)
                                            .attr("cy", y2)
                                            .attr("class", "arrowlineEnd")
                                            .attr("opacity", 0)
                                            .call(d3.drag()
                                            .on("drag", function(d){
                                                var obj = d3.select(this);
                                                obj.style("cursor", "move");
                                                var coords = d3.mouse(svg.node());
                                                var posx = Math.round(xscale().invert(coords[0] - margin.left));
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
                                                          .attr("x2", xscale()(posx)+margin.left)
                                                          .attr("y2", posy)
                                                          .attr("ref", "arrow"+posx+"_"+m.ref[2]+"_"+posy+"_"+m.ref[3]);
                                                    obj.attr("cx", xscale()(posx)+margin.left)
                                                       .attr("cy", posy);
                                                }
                                            })
                                            .on("end", function(d){
                                                 d3.select(this).style("cursor", "default");
                                                 var ref=d3.select(this.parentNode).attr("ref");
                                                 addNewHistory({
													undo:function(){
														delete(x.markers[ref]);
														var m={
															"ref" : [oldstart.x1, oldstart.y1, oldstart.x2, oldstart.y2],
															"markertype" : 3,
															"color" : oldstart.col,
															"opacity" : 1,
															"linetype" : "solid",
															"linewidth" : 1
														};
														x.markers["arrow"+oldstart.x1+"_"+oldstart.x2+"_"+oldstart.y1+"_"+oldstart.y2]=m;
														self.redraw();
													},
													redo:function(){}
												});
                                            })
                                            .on("start", function(d){
                                            	d3.select(this).style("cursor", "move");
                                            	oldstart={x1:x1, y1:y1, x2:x2, y2:y2, col:m.color};
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
                                     .on("mouseover", function(){ d3.select(this).attr("fill", "#DEDEDE").attr("opacity", 0.5); })
                                     .on("mouseout", function(){ d3.select(this).attr("fill", "white").attr("opacity", 0); })
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
                            l= self.g.append("rect")
									   .attr("stroke", "none")
									   .attr("fill", "white")
									   .attr("opacity", 0)
									   .attr("x", xscale()(m.ref) + margin.left - 5)
									   .attr("width", 10)
									   .attr("y", margin.top)
									   .attr("height", +svg.attr("height") - margin.bottom - margin.top)
                                       .style("cursor", "move")
                                       .attr("ref", "line"+m.ref)
                                       .attr("class", "Marker");
                            var ll=self.g.append("line")
                                       .attr("y1", margin.top)
                                       .attr("y2", +svg.attr("height") - margin.bottom)
                                       .attr("x1", xscale()(m.ref) + margin.left)
                                       .attr("x2", xscale()(m.ref) + margin.left)
                                       .attr("stroke", m.color)
                                       .attr("stroke-width", "1px")
                                       .attr("fill", "none")
                                       .style("cursor", "move")
                                       .attr("ref", "line"+m.ref)
                                       .attr("class", "Marker");
                            if(m.linetype==="dashed"){
                                ll.style("stroke-dasharray", ("3, 3"));
                            }else{
                                //solid
                            }
                            l.call(d3.drag().on("drag", self.draggedLine)
                                            .on("end", self.dragendLine));
                            break;
                        case 1:
                            var newWid=xscale()(m.ref+m.linewidth) - xscale()(m.ref);
                            if(newWid>1){
                              l= self.g.append("rect")
                                   .attr("id", "rect"+m.ref)
                                   .attr("stroke", "none")
                                   .attr("fill", m.color)
                                   .attr("x", xscale()(m.ref) + margin.left)
                                   .attr("width", newWid)
                                   .attr("y", margin.top)
                                   .attr("height", +svg.attr("height") - margin.bottom - margin.top)
                                   .attr("ref", "rect"+m.ref)
                                   .style("opacity", m.opacity)
                                   .attr("class", "Marker");
                            self.g.append("line")
                                   .attr("id", "rectL"+m.ref)
                                   .attr("y1", margin.top)
                                   .attr("y2", +svg.attr("height") - margin.bottom)
                                   .attr("x1", xscale()(m.ref) + margin.left)
                                   .attr("x2", xscale()(m.ref) + margin.left)
                                   .attr("stroke", "white")
                                   .attr("stroke-width", "2px")
                                   .style("opacity", 0)
                                   .style("cursor", "ew-resize")
                                   .attr("ref", "rect"+m.ref)
                                   .call(d3.drag().on("drag", self.resizeRectL).on("end", self.resizeRectLEnd));
                            self.g.append("line")
                                   .attr("id", "rectR"+m.ref)
                                   .attr("y1", margin.top)
                                   .attr("y2", +svg.attr("height") - margin.bottom)
                                   .attr("x1", xscale()(m.ref+m.linewidth) + margin.left)
                                   .attr("x2", xscale()(m.ref+m.linewidth) + margin.left)
                                   .attr("stroke", "white")
                                   .attr("stroke-width", "2px")
                                   .style("opacity", 0)
                                   .style("cursor", "ew-resize")
                                   .attr("ref", "rect"+m.ref)
                                   .call(d3.drag().on("drag", self.resizeRectR).on("end", self.resizeRectREnd));
                            }
                            break;
                        case 2:
							var opt = textDefaultOptions();
							opt.id = "text"+m.ref[0] + "_" + m.ref[1];
							opt.ref = "text"+m.ref[0] + "_" + m.ref[1];
							opt.color = m.color;
							opt.x = xscale()(m.ref[0]) + margin.left;
							opt.y = m.ref[1];
							opt.anchor = "start";
							opt.fontsize = m.fontsize;
							opt.angle = m.angle;
							opt.text = m.text;
							opt.vp = self.g;
							opt.cls = "Mlabel";
							opt.colorPickerId = 9;
							l = new Label(opt);
                            break;
                        case 3:
                            l= new Arrow(self.g, m);
                            break;
                        case 4:
                        	m.vp = self.g;
                        	var llab = new legendLabel(m);
                        	break;
                    }
                    if(m.markertype<2){
                        l.on("click", function(){
                        	var obj = d3.select(this);
							var picked = function(col){
								var tmpstatus = {k:obj.attr("ref"), v:x.markers[obj.attr("ref")].color};
								addNewHistory({
									undo:function(){
										x.markers[tmpstatus.k].color = tmpstatus.v;
										plotregion.renew();
									},
									redo:function(){}
								});
								x.markers[obj.attr("ref")].color = col;
								self.redraw();
							};
							ColorPicker(this, picked);
                        })
                         .on("dblclick", function(){
                            var obj = d3.select(this);
                         	var tmpstatus = {
                         		k:obj.attr("ref"),
                         		v:clone(x.markers[obj.attr("ref")])
                         	};
                         	addNewHistory({
                         		undo:function(){
                         			x.markers[tmpstatus.k] = clone(tmpstatus.v);
                         			self.redraw();
                         		},
                         		redo:function(){}
                         	});
                            delete x.markers[obj.attr("ref")];
                            obj.remove();
                            self.redraw();
                            self.svgDefault = false;
                         });
                     }
                }
            }
            self.redraw();
            return(self);
        };
        
        // Ruler
        var currentLayer = 0;
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
            /*self.yrule_label2 = svg.insert("text", ":first-child")
                                    .attr("fill", "#333")
                                    .attr("y", 10)
                                    .attr("x", +svg.attr("width")-5)
                                    .attr("text-anchor", "end")
                                    .style("font-size", "0.75em")
                                    .text("");*/
            self.trackNum = function(xpos, ypos){
                var H = (ypos - margin.top - 40)/heightF();
                if(H>1) return(ypos);
                for(var i=0; i<trackNames().length; i++){
                    H = H - x.height[trackNames()[i]];
                    if(H<0){
						currentLayer = i;
                        if(x.type[trackNames()[i]]==="data"){
                            var a=x.tracklist[trackNames()[i]].dat,
                                b=x.tracklist[trackNames()[i]].dat2;
                            var c="";
                            if(a.length) c = d3.format(".2f")(a[xpos]);
                            if(b.length) c = c + "; "+d3.format(".2f")(b[xpos]);
                            return(c);
                        }else{
                            return(trackNames()[i]);
                        }
                    }
                }
				currentLayer = i;
                return(trackNames()[i]);
            }
            self.rulemove = function(coords){
                var xpos = Math.round(xscale().invert(coords[0]-margin.left));
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
                //self.yrule_label2.attr("y", coords[1]-5)
                //           .text(ypos);
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
        };
        
        //ColorPicker
        var cp;//colorPicker;
        var color = "#000";
        var cpCheckAll = false;
        var ColorPicker = function (target, picked, allbutton=false) {
        	var self = this;
        	var target = d3.select(target);
        	var currentId = Number(target.attr("kvalue"));
            var colorScale = ["#FFD300","#FFFF00","#A2F300","#00DB00","#00CD00","#00FFFF",
                              "#00B7FF","#0000FF","#1449C4","#4117C7","#820AC3","#DB007C",
                              "#FF0000","#FF00FF","#FF7400","#FFAA00"];//change to color blindness safe?
            var getColor = function (i) {
                return colorScale[i];
            };
            if(typeof(target.attr("fill"))!="undefined"){
            	color = target.attr("fill");
            }else{
				if(typeof(target.attr("stroke"))!="undefined"){
					color = target.attr("stroke");
				}
            }
            defaultColor = color || getColor(0);
        
            self.pickedColor = defaultColor;
            self.defaultPicked = function (col) {
            	switch(picked){
            		case 0: //data track dat
            			x.tracklist[trackNames()[currentId]].style.color[0] = col;
            			break;
            		case 1: //data track dat2
            			x.tracklist[trackNames()[currentId]].style.color[1] = col;
            			break;
            		case 2: //track baseline
            			//x.tracklist[trackNames()[currentId]].style.color = col;
            			var featureLayerID = Number(target.attr("ref"));
            			for(var i=0; i<x.tracklist[trackNames()[currentId]].dat.featureLayerID.length; i++){
            				if(x.tracklist[trackNames()[currentId]].dat.featureLayerID[i]==featureLayerID){
            					x.tracklist[trackNames()[currentId]].dat.fill[i] = col;
            				}
            			}
            			break;
            		case 3: //label or lines
            			var tmpstatus = {k:currentId,v:x.color[trackNames()[currentId]]};
            			addNewHistory({
            				undo:function(){
            					x.color[trackNames()[tmpstatus.k]] = tmpstatus.v;
            					plotregion.renew();
            				},
            				redo:function(){}
            			});
            			x.color[trackNames()[currentId]] = col;
            			break;
            		case 4: //lollipop nodes
            			var poskey = Number(target.attr("poskey"));
            			var datatrack = target.attr("datatrack");
            			var k = Number(target.attr("kvalue"));
            			x.tracklist[trackNames()[k]][datatrack].color[poskey]=col;
            			break;
            		case 5: // lollipop lines
            			var poskey = Number(target.attr("poskey"));
            			var datatrack = target.attr("datatrack");
            			var k = Number(target.attr("kvalue"));
            			x.tracklist[trackNames()[k]][datatrack].border[poskey]=col;
            			break;
            		case 6: // lollipop label
            			var poskey = Number(target.attr("poskey"));
            			var datatrack = target.attr("datatrack");
            			var k = Number(target.attr("kvalue"));
            			var tmpstatus = {k:k,datatrack:datatrack,poskey:poskey,v:x.color["lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+poskey]};
            			addNewHistory({
            				undo:function(){
            					x.color["lolliplotTrackLabel_"+tmpstatus.k+"_"+tmpstatus.datatrack+"_"+tmpstatus.poskey] = tmpstatus.v;
            					plotregion.renew();
            				},
            				redo:function(){}
            			});
            			x.color["lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+poskey] = col;
            			break;
            		case 7: // lollipop baseline bottom line
            			break;
            		case 8: // lollipop baseline top line
            			break;
            		case 9: //markers
            			var tmpstatus = {k:target.attr("ref"),v:x.markers[target.attr("ref")].color};
            			addNewHistory({
            				undo:function(){
            					x.markers[tmpstatus.k].color = tmpstatus.v;
            					plotregion.renew();
            				},
            				redo:function(){}
            			});
            			x.markers[target.attr("ref")].color = col;
            			break;
            		default: 
            			console.log(target);
            	}
            	
                color = col;
                plotregion.renew();
                newCP();//keep it on top
            };
            var clicked = function () {
            	if(typeof(picked)=="function"){
            		picked(self.pickedColor);
            	}else{
                	self.defaultPicked(self.pickedColor);
                }
            };
        
            var pie = d3.pie().sort(null);
            var arc = d3.arc().innerRadius(25).outerRadius(50);
            var currentCoor = [d3.event.x, d3.event.y];
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
						.attr("fill", "#000")
						.attr("stroke", "#fff")
						.attr("stroke-width", 2)
						.attr("r", 10)
						.attr("cx", -45)
						.attr("cy", 45)
						.on("mouseover", function () {
							var fill = target.attr("fill");
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
			
				var colLen = [];
				for(var i=0; i<colorScale.length; i++){
					colLen.push(1);
				}
				cp.datum(colLen)
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
				var frm = cp.append("foreignObject")
							.attr("x", -28)
							.attr("y", 50)
							.attr("width", 50)
							.attr("height", 20);
				var inp = frm.append("xhtml:form")
								.append("input")
								.attr("value", defaultColor)
								.attr("style", "width:50px;")
								.on("keypress", function(){
									// IE fix
									if (!d3.event)
										d3.event = window.event;
									var e = d3.event;
									if (e.keyCode == 13)
									{
										if (typeof(e.cancelBubble) !== 'undefined') // IE
										  e.cancelBubble = true;
										if (e.stopPropagation)
										  e.stopPropagation();
										e.preventDefault();
										 var fill = inp.node().value;
										 if(/^#(?:[0-9a-fA-F]{3}){1,2}$/.exec(fill)){
											 self.pickedColor = fill;
											 plate.attr("fill", self.pickedColor);
											 clicked();
										 }
									}
								});
				if(allbutton){
					var check4all = cp.append("g")
										.attr("width", 50)
										.attr("height", 20)
										.attr("transform", "translate(-45 -45)");
					var frm2 = check4all.append("foreignObject")
								.attr("x", -14)
								.attr("y", -10)
								.attr("width", 20)
								.attr("height", 20);
					var inp2 = frm2.append("xhtml:form")
									.append("input")
									.attr("type", "checkbox")
									.attr("name", "check4all")
									.attr("value", "check4all")
									.on("change", function(){
										cpCheckAll = d3.select(this).property("checked");
									});
					if(cpCheckAll) inp2.attr("checked", true);
					var check4allLabel = check4all.append("text")
												  .attr("fill", "#000")
												  .attr("x", 3)
												  .attr("y", 5)
												  .text("all")
												  .style("cursor", "default");
				}
            };
            newCP();
            return(self);
        };
        
        // context menu
        function contextMenu(mg){
			var height = 24,
				width = 150, 
				margin = 10, // fraction of width
				x0 = 0,
				y0 = 0,
				items = [
					{label:'add arrow label',
					 onMouseClick: function(x0, y0){mg.addArrowLabel(x0, y0);},
					 check: function(){return true;}
					}, 
					{
					 label:'add legend label',
					 onMouseClick: function(x0, y0) {mg.addLegendLabel(x0, y0)},
					 check: function(){return true;}
					},
					{
					 label:'remove all guidelines',
					 onMouseClick: function(){
					 	var tmpstatus = {};
					 	var keys = Object.keys(x.markers);
					 	for(var i=0; i<keys.length; i++){
					 		if(x.markers[keys[i]].markertype==0){
					 			d3.select('.Marker[ref="'+keys[i]+'"]').remove();
					 			tmpstatus[keys[i]] = clone(x.markers[keys[i]]);
					 			delete x.markers[keys[i]];
					 		}
					 	}
					 	plotregion.renew();
					 	addNewHistory({
										undo:function(){
											var keys = Object.keys(tmpstatus);
											for(var i=0; i<keys.length; i++){
												x.markers[keys[i]] = clone(tmpstatus[keys[i]]);
											}
											plotregion.renew();
										},
										redo:function(){}
									  });
					 },
					 check: function(){return true;}
					},
					{label: 'decrease gene track height',
					  check: function(){
					  	var tx=false;
					  	var dat=false;
					  	var keys=Object.keys(x.type);
					  	for(var i=0; i<keys.length; i++){
					  		if(x.type[keys[i]]=="data" || x.type[keys[i]]=="lollipopData"){
					  			dat=true;
					  		}else{
					  			tx=true;
					  		}
					  	}
					  	return(dat && tx);
					  },
					  onMouseClick: function(){
					  	var tmpstatus = clone(x.height);
					  	addNewHistory({
										undo:function(){
											x.height = tmpstatus;
											plotregion.renew();
										},
										redo:function(){}
									  });
					  	d3.selectAll('g[type = "transcript"]')
							 .each(function(){
								var obj=d3.select(this);
								x.height[trackNames()[obj.attr('kvalue')]] = 
								  x.height[trackNames()[obj.attr('kvalue')]]*.9;
							 });
					  	d3.selectAll('g[type = "gene"]')
							 .each(function(){
								var obj=d3.select(this);
								x.height[trackNames()[obj.attr('kvalue')]] = 
								  x.height[trackNames()[obj.attr('kvalue')]]*.9;
							 });
					  	var totalH=0;
					  	for(var i=0; i<trackNames().length;i++){
					  		totalH+=x.height[trackNames()[i]];
					  	}
					  	for(var i=0; i<trackNames().length; i++){
						   x.height[trackNames()[i]] = x.height[trackNames()[i]]/totalH;
					    }
					  	plotregion.renew();
					  }
					}, 
					{label: 'increase gene track height',
					  check: function(){
					  	var tx=false;
					  	var dat=false;
					  	var keys=Object.keys(x.type);
					  	for(var i=0; i<keys.length; i++){
					  		if(x.type[keys[i]]=="data" || x.type[keys[i]]=="lollipopData"){
					  			dat=true;
					  		}else{
					  			tx=true;
					  		}
					  	}
					  	return(dat && tx);
					  },
					  onMouseClick: function(){
					  	var tmpstatus = clone(x.height);
					  	addNewHistory({
										undo:function(){
											x.height = tmpstatus;
											plotregion.renew();
										},
										redo:function(){}
									  });
					  	d3.selectAll('g[type = "transcript"]')
							 .each(function(){
								var obj=d3.select(this);
								x.height[trackNames()[obj.attr('kvalue')]] = 
								  x.height[trackNames()[obj.attr('kvalue')]]*1.1;
							 });
					  	d3.selectAll('g[type = "gene"]')
							 .each(function(){
								var obj=d3.select(this);
								x.height[trackNames()[obj.attr('kvalue')]] = 
								  x.height[trackNames()[obj.attr('kvalue')]]*1.1;
							 });
					  	var totalH=0;
					  	for(var i=0; i<trackNames().length;i++){
					  		totalH+=x.height[trackNames()[i]];
					  	}
					  	for(var i=0; i<trackNames().length; i++){
						   x.height[trackNames()[i]] = x.height[trackNames()[i]]/totalH;
					    }
					  	plotregion.renew();
					  }
					},
					{label: 'change Y-axis position',
					 check:function(){
					 	if(eventLayer=="") return false;
					  	return(x.type[eventLayer]=="data");
					  },
					  onMouseClick: function(){
					  	var tmpstatus = {k:eventLayer,
					  					 v:x.tracklist[eventLayer].style.yaxis.main};
					  	addNewHistory({
										undo:function(){
											x.tracklist[tmpstatus.k].style.yaxis.main = tmpstatus.v;
											plotregion.renew();
										},
										redo:function(){}
									  });
					  	x.tracklist[eventLayer].style.yaxis.main = !x.tracklist[eventLayer].style.yaxis.main;
					  	plotregion.renew();
					  }
					},
					{label: 'condense lollipops',
					 check: function(){
					 	if(eventLayer=="") return false;
					 	if(typeof(x.tracklist[eventLayer].dat2.type)!="undefined"){
							if(x.tracklist[eventLayer].dat2.type[0] == "dandelion"){
								return false;
							}
						}
						return(x.type[eventLayer]!="data");
					 },
					 onMouseClick: function(){
						if(x.type[eventLayer]!="data"){
							if(x.type[eventLayer]=="lollipopData"){
								condenseLollipops(x.tracklist[eventLayer].dat, true);
							}
							if(typeof(x.tracklist[eventLayer].dat2.labpos)!="undefined"){
								condenseLollipops(x.tracklist[eventLayer].dat2, true);
							}
							plotregion.renew();
					 	}
					 }
					},
					{label: 'expand lollipops',
					 check: function(){
					 	if(eventLayer=="") return false;
					 	if(typeof(x.tracklist[eventLayer].dat2.type)!="undefined"){
							if(x.tracklist[eventLayer].dat2.type[0] == "dandelion"){
								return false;
							}
						}
						return(x.type[eventLayer]!="data");
					 },
					 onMouseClick: function(){
						if(x.type[eventLayer]!="data"){
							if(x.type[eventLayer]=="lollipopData"){
								condenseLollipops(x.tracklist[eventLayer].dat, false);
							}
							if(typeof(x.tracklist[eventLayer].dat2.labpos)!="undefined"){
								condenseLollipops(x.tracklist[eventLayer].dat2, false);
							}
							plotregion.renew();
					 	}
					 }
					},
					{label: 'rotate all labels 45º(CW)',
					 check: function(){
					 	if(eventLayer=="") return false;
					 	if(x.type[eventLayer]!="data"){
					 		if(x.type[eventLayer]=="lollipopData"){
					 			return(typeof(x.tracklist[eventLayer].dat["textlabel"])!="undefined");
					 		}else{
					 			return(typeof(x.tracklist[eventLayer].dat2["textlabel"])!="undefined");
					 		}
					 	}else{
					 		return false;
					 	}
					 },
					 onMouseClick: function(){
					 	var tmpstatus = {k:eventLayer};
					 	if(x.type[eventLayer]=="lollipopData"){
							tmpstatus = {k:eventLayer,
										 d:"dat",
										 v:clone(x.tracklist[eventLayer].dat["label.parameter.rot"])};
							if(typeof(x.tracklist[eventLayer].dat["label.parameter.rot"])=="undefined"){
								x.tracklist[eventLayer].dat["label.parameter.rot"] = [];
							}
					 		for(var i=0; i<x.tracklist[eventLayer].dat["textlabel"].length; i++){
					 			if(typeof(x.tracklist[eventLayer].dat["label.parameter.rot"][i])=="undefined"){
					 				x.tracklist[eventLayer].dat["label.parameter.rot"][i] = 90;
					 			}
								x.tracklist[eventLayer].dat["label.parameter.rot"][i]-=45;
							}
						}else{
							tmpstatus = {k:eventLayer,
										 d:"dat2",
										 v:clone(x.tracklist[eventLayer].dat2["label.parameter.rot"])};
							if(typeof(x.tracklist[eventLayer].dat2["label.parameter.rot"])=="undefined"){
								x.tracklist[eventLayer].dat["label.parameter.rot"] = [];
							}
							for(var i=0; i<x.tracklist[eventLayer].dat2["textlabel"].length; i++){
					 			if(typeof(x.tracklist[eventLayer].dat2["label.parameter.rot"][i])=="undefined"){
					 				x.tracklist[eventLayer].dat2["label.parameter.rot"][i] = 90;
					 			}
								x.tracklist[eventLayer].dat2["label.parameter.rot"][i]-=45;
							}
						}
						plotregion.renew();
					 	addNewHistory({
					 		undo:function(){
					 			x.tracklist[tmpstatus.k][tmpstatus.d]["label.parameter.rot"] = clone(tmpstatus.v);
					 			plotregion.renew();
					 		},
					 		redo:function(){}
					 	});
					 }
					},
					{label: 'change lollipops type',
					 check: function(){
					 	if(eventLayer=="") return false;
					 	if(x.type[eventLayer]!="data"){
					 		if(x.type[eventLayer]=="lollipopData"){
					 			return(typeof(x.tracklist[eventLayer].dat["stack.factor"])=="undefined");
					 		}else{
					 			return(typeof(x.tracklist[eventLayer].dat2["stack.factor"])=="undefined");
					 		}
					 	}else{
					 		return false;
					 	}
					 },
					 onMouseClick: function(){
					 	var types = {"circle":"pin", "pin":"pie", "pie":"dandelion", "dandelion":"circle"};
					 	var tmpstatus = {k:eventLayer};
						if(x.type[eventLayer]!="data"){
							if(x.type[eventLayer]=="lollipopData"){
								tmpstatus = {k:eventLayer,
												 d:"dat",
												 v:clone(x.tracklist[eventLayer].dat.type)};
								if(typeof(x.tracklist[eventLayer].dat["stack.factor"])=="undefined"){
									if(typeof(x.tracklist[eventLayer].dat.type)!="undefined"){
										for(var i=0; i<x.tracklist[eventLayer].dat.type.length; i++){
											x.tracklist[eventLayer].dat.type[i] = 
												types[x.tracklist[eventLayer].dat.type[i]];
										}
									}else{
										x.tracklist[eventLayer].dat.type = [];
										for(var i=0; i<x.tracklist[eventLayer].dat.labpos.length; i++){
											x.tracklist[eventLayer].dat.type.push("pin");
										}
									}
								}
							}else{
								tmpstatus = {k:eventLayer,
												 d:"dat2",
												 v:clone(x.tracklist[eventLayer].dat2.type)};
								if(typeof(x.tracklist[eventLayer].dat2["stack.factor"])=="undefined"){
									if(typeof(x.tracklist[eventLayer].dat2.labpos)!="undefined"){
										if(typeof(x.tracklist[eventLayer].dat2.type)!="undefined"){
											for(var i=0; i<x.tracklist[eventLayer].dat2.type.length; i++){
												x.tracklist[eventLayer].dat2.type[i] = 
													types[x.tracklist[eventLayer].dat2.type[i]];
											}
										}else{
											x.tracklist[eventLayer].dat2.type = [];
											for(var i=0; i<x.tracklist[eventLayer].dat2.labpos.length; i++){
												x.tracklist[eventLayer].dat2.type.push("pin");
											}
										}
									}
								}
							}
						}
						plotregion.renew();
					 	addNewHistory({
					 		undo:function(){
					 			x.tracklist[tmpstatus.k][tmpstatus.d].type = clone(tmpstatus.v);
					 			plotregion.renew();
					 		},
					 		redo:function(){}
					 	});
					 }
					},
					{label:'condense isoforms', 
					check: function(){
					 	if(eventLayer=="") return false;
						return(x.type[eventLayer]=="gene" || x.type[eventLayer]=="transcript")
					},
					onMouseClick: function(){
						var tmpstatus = {k:eventLayer,v:parameter.isoformR[eventLayer]};
						addNewHistory({
					 		undo:function(){
					 			parameter.isoformR[tmpstatus.k] = tmpstatus.v;
					 			plotregion.renew();
					 		},
					 		redo:function(){}
					 	});
						parameter.isoformR[eventLayer] +=  0.25;
						plotregion.renew();
					}},
					{label:'expand isoforms', 
					check: function(){
					 	if(eventLayer=="") return false;
						return(x.type[eventLayer]=="gene" || x.type[eventLayer]=="transcript")
					},
					onMouseClick: function(){
						var tmpstatus = {k:eventLayer,v:parameter.isoformR[eventLayer]};
						addNewHistory({
					 		undo:function(){
					 			parameter.isoformR[tmpstatus.k] = tmpstatus.v;
					 			plotregion.renew();
					 		},
					 		redo:function(){}
					 	});
						parameter.isoformR[eventLayer] -= 0.25;
						plotregion.renew();
					}},
					{label:'increase maxgap', 
					check: function(){
					 	if(eventLayer=="") return false;
						if(typeof(x.tracklist[eventLayer].dat2.type)!="undefined"){
							if(x.tracklist[eventLayer].dat2.type[0] == "dandelion"){
								return true;
							}
						}
						return false;
					},
					onMouseClick: function(){
						var tmpstatus = parameter.maxGAP;
						addNewHistory({
					 		undo:function(){
					 			parameter.maxGAP = tmpstatus;
					 			plotregion.renew();
					 		},
					 		redo:function(){}
					 	});
						parameter.maxGAP *= 1.25;
						plotregion.renew();
					}},
					{label:'decrease maxgap', 
					check: function(){
					 	if(eventLayer=="") return false;
						if(typeof(x.tracklist[eventLayer].dat2.type)!="undefined"){
							if(x.tracklist[eventLayer].dat2.type[0] == "dandelion"){
								return true;
							}
						}
						return false;
					},
					onMouseClick: function(){
						var tmpstatus = parameter.maxGAP;
						addNewHistory({
					 		undo:function(){
					 			parameter.maxGAP = tmpstatus;
					 			plotregion.renew();
					 		},
					 		redo:function(){}
					 	});
						parameter.maxGAP *= .75;
						plotregion.renew();
					}}/*,
					{label:'flip coordinates', 
					check: function(){},
					onMouseClick: function(){
						alert("coming soon.");
					}}*/
				];
	        var eventLayer = trackNames()[currentLayer];
	        var condenseLollipops = function(data, condense=true){
	        	var labpos = [];
	        	var labCenter = 0;
	        	var tmpstatus = clone(data.labpos);
	        	addNewHistory({
	        		undo:function(){
	        			data.labpos = clone(tmpstatus);
	        			plotregion.renew();
	        		},
	        		redo:function(){}
	        	});
	        	for(var i=0; i<data.labpos.length; i++){
	        		labpos.push(data.labpos[i]);
	        		labCenter +=data.labpos[i];
	        	}
	        	labCenter = labCenter/labpos.length;
	        	var labMin = Math.min(...labpos);
	        	var labMax = Math.max(...labpos);
	        	var labShift = (labMax - labMin)/(labpos.length+1)/10;
	        	var labCp = 0;
	        	var labCd = labMax;
	        	for(var i=0; i<data.labpos.length; i++){
	        		var thisLabCd = Math.abs(data.labpos[i] - labCenter);
	        		if(thisLabCd < labCd){
	        			labCd = thisLabCd;
	        			labCp = i;
	        		}
	        	}
	        	if(typeof(data["stack.factor"])=="undefined"){
					for(var i=0; i<data.labpos.length; i++){
						if(condense){
							data.labpos[i] += (labCp-i)*labShift;
						}else{
							data.labpos[i] -= (labCp-i)*labShift;
						}
					}
	        	}else{
	        		//group by labpos
	        		var labposGroup = {};
	        		var labposGroupInv = [];
	        		for(var i=0; i<data.labpos.length; i++){
	        			if(typeof(labposGroup["X"+data.labpos[i]])=="undefined"){
	        				labposGroup["X"+data.labpos[i]] = [];
	        			}
						labposGroup["X"+data.labpos[i]].push(i);
						labposGroupInv.push("X"+data.labpos[i]);
					}
					var key = Object.keys(labposGroup);
	        		for(var i=0; i<key.length; i++){
						for(var j=0; j<labposGroup[key[i]].length; j++){
							if(condense){
								data.labpos[labposGroup[key[i]][j]] += (key.indexOf(labposGroupInv[labCp])-i)*labShift;
							}else{
								data.labpos[labposGroup[key[i]][j]] -= (key.indexOf(labposGroupInv[labCp])-i)*labShift;
							}
						}
	        		}
	        	}
	        };
			function menu(x, y) {
				x0=x;
				y0=y;
				d3.select('.context-menu').remove();
			    eventLayer = trackNames()[currentLayer];
			    if(typeof(eventLayer)=="undefined"){
			    	eventLayer = "";
			    }
			    if(y + height*items.length > +svg.attr("height")){
			    	y = +svg.attr("height") - height*items.length;
			    }
			    if(x + width > +svg.attr("width")){
			    	x = +svg.attr("width") - width;
			    }
				// Draw the menu
				svg.append('g').attr('class', 'context-menu')
					.selectAll('.menu-entry')
					.data(items.filter(function(d){return d.check();})).enter()
					.append('g').attr('class', 'menu-entry')
					.style('cursor', 'pointer')
					.on('mouseover', function(){ 
						d3.select(this).select('rect')
						  .style('fill', '#EEEEEE') 
					})
					.on('mouseout', function(){ 
						d3.select(this).select('rect')
						.style('fill', '#FAFAFA')
						.style('stroke', 'white')
						.style('stroke-width', '1px') 
					});
		
				d3.selectAll('.menu-entry')
					.append('rect')
					.attr('x', x)
					.attr('y', function(d, i){ return y + (i * height); })
					.attr('width', width)
					.attr('height', height)
					.style('fill', '#FAFAFA')
					.style('stroke', 'white')
					.style('stroke-width', '1px')
					.on('click', function(d){
						d.onMouseClick();
					});
		
				d3.selectAll('.menu-entry')
					.append('text')
					.text(function(d){ return d.label; })
					.attr('x', x)
					.attr('y', function(d, i){ return y + (i * height); })
					.attr('dy', height - margin / 2)
					.attr('dx', margin)
					.style('fill', function(d, i){ 
						return d.check()?'darkblue':'gray';
					})
					.style('font-size', '13')
					.on('click', function(d){
						d.onMouseClick(x0, y0);
					});

				// Other interactions
				d3.select('body')
					.on('click', function() {
						d3.select('.context-menu').remove();
					});

			}

			return menu;
        }
        
        // make the plot region resizable
        var Margin = function(plotregion){
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
								.on("mouseover", function(){d3.select(this).attr("stroke", 'gray').style("opacity", 0.5)})
								.on("mouseout", function(){d3.select(this).attr("stroke", 'white').style("opacity", 0)})
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
								.on("mouseover", function(){d3.select(this).attr("stroke", 'gray').style("opacity", 0.5)})
								.on("mouseout", function(){d3.select(this).attr("stroke", 'white').style("opacity", 0)})
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
								.on("mouseover", function(){d3.select(this).attr("stroke", 'gray').style("opacity", 0.5)})
								.on("mouseout", function(){d3.select(this).attr("stroke", 'white').style("opacity", 0)})
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
								.on("mouseover", function(){d3.select(this).attr("stroke", 'gray').style("opacity", 0.5)})
								.on("mouseout", function(){d3.select(this).attr("stroke", 'white').style("opacity", 0)})
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
                    case 1://marginbottom
                        dy = +svg.attr("height") - dy;
                        if(dy<0) dy=0;
                        margin.bottom = dy;
                        plotregion.g.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
                        plotregion.renew();
                        break;
                    case 2://marginleft
                        if(dx<0) dx=0;
                        margin.left = dx;
                        plotregion.g.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
                        plotregion.renew();
                        break;
                    case 3://margintop
                        if(dy<0) dy=0;
                        margin.top = dy;
                        plotregion.g.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
                        plotregion.renew();
                        break;
                    case 4://marginright
                        dx = +svg.attr("width") - dx;
                        if(dx<0) dx=0;
                        margin.right = dx;
                        plotregion.renew();
                        break;
                }
            }
            return(self);
        };
        
        
        //arrow of intron
        parameter.arrowId = 0;
        var plotArrow = function(container, data, start, end, ypos, strand, color, k){
            if(strand==="*") return(null);
            parameter.arrowId +=1;
            var marker = container.append("marker")
                     .attr("id", "arrow"+safeNames()[k]+"_"+parameter.arrowId)
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
                     .attr("orient", strand==="-"?180:0)
                     .append("path")
                     .attr("d", "M 0 0 L 10 5 L 0 10");
            for(var i=0; i<data.length; i++){
                for(var j=data[i].x0; j<=data[i].x1; j+=10){
                    if(j > start && j < end){
                        container.append("line")
                            .attr("x1",  j)
                            .attr("y1", ypos)
                            .attr("x2",  j)
                            .attr("y2", ypos)
                            .attr("marker-end", "url(#arrow"+safeNames()[k]+"_"+parameter.arrowId+")")
                            .attr("stroke-width", 1)
                            .attr("class", "geneArrow track"+safeNames()[k])
                            .attr("kvalue", k);
                    }
                }
            }
            
        }
        
        // lolliplot
        parameter.circleR=[];
        parameter.YposSTART = [];
        parameter.maxGAP=1;
        var lolliplot = function(lolli, trackdat, xscale, yscale, k, datatrack, ypos, from){
        	var pos;
        	var defaultCircleR=.05;
        	if(typeof(parameter.circleR[safeNames()[k]+"_"+datatrack])=="undefined"){
        		parameter.circleR[safeNames()[k]+"_"+datatrack] = defaultCircleR;
        	}
        	var thisCR=parameter.circleR[safeNames()[k]+"_"+datatrack];
        	var thisCR2 = 1 - 2*thisCR;
        	thisCR = 1 - thisCR;
        	var yposStart = 0;
        	if(typeof(parameter.YposSTART[safeNames()[k]+"_"+datatrack])!="undefined") yposStart = parameter.YposSTART[safeNames()[k]+"_"+datatrack];
        	switch(ypos){
        		case .165:
        			pos=[.165, .35+yposStart, .45+yposStart*2, .55+yposStart*3, .6+yposStart*3, thisCR, thisCR2];
        			break;
        		case .45:
        			pos=[.45, .55+yposStart, .65+yposStart*2, .70+yposStart*3, .75+yposStart*3, thisCR, thisCR2];
        			break;
        		case .55:
        			pos=[.55, .45-yposStart, .35-yposStart*2, .30-yposStart*3, .25-yposStart*3, thisCR, thisCR2];
        			break;
        		case .825:
        			pos=[.825, .65-yposStart, .55-yposStart*2, .50-yposStart*3, .45-yposStart*3, thisCR, thisCR2];
        			break;
        	}
        	function dragstarted(d) {
			  d3.select(this).style("cursor", "move").raise().classed("active", true);
			  var ele = d3.select(this);
			  var k = Number(ele.attr("kvalue"));
			  var datatrack = ele.attr("datatrack");
			  var poskey = Number(ele.attr("poskey"));
			  var oldpos = x.tracklist[trackNames()[k]][datatrack].labpos[poskey];
			  addNewHistory({
			  	undo:function(){
			  		x.tracklist[trackNames()[k]][datatrack].labpos[poskey]=oldpos;
			  		plotregion.renew();
			  	},
			  	redo:function(){}
			  });
			}
			function dragended(d) {
			  d3.select(this).style("cursor", "default").classed("active", false);
			}
			function draggedGroup(d){
			  var ele = d3.select(this);
			  var k = Number(ele.attr("kvalue"));
			  var datatrack = ele.attr("datatrack");
			  var poskey = Number(ele.attr("poskey"));
			  var posx = xscale.invert(d3.mouse(svg.node())[0]-margin.left);
			  x.tracklist[trackNames()[k]][datatrack].labpos[poskey]=posx;
			  //self renew;
			  d3.select("#lollipopNodeLinker0_"+safeNames()[k]+"_"+datatrack+"_"+poskey)
			  	.attr("x2", xscale(posx));
			  d3.select("#lollipopNodeLinker1_"+safeNames()[k]+"_"+datatrack+"_"+poskey)
			  	.attr("x1", xscale(posx))
			  	.attr("x2", xscale(posx));
			  
			  var type=x.type[trackNames()[k]];
			  if(typeof(x.tracklist[trackNames()[k]][datatrack]["stack.factor"])!="undefined"){
					type = "pie.stack";
				}
			  if(typeof(x.tracklist[trackNames()[k]][datatrack].type)!="undefined"){
			  	type = x.tracklist[trackNames()[k]][datatrack].type[poskey];
			  }
			  switch(type){
			  	case "gene":
			  	case "transcript":
			  	case "lollipopData":
			  	case "circle":
			  	case "pin":
				  d3.select("#lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+poskey)
					.selectAll("circle").attr("cx", xscale(posx));
				  d3.select("#lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+poskey)
					.select("path").attr("d", "m "+(xscale(posx)-2*yscale(pos[5]))+" "+(circenter+((ypos<.5?-1:1)*2*yscale(pos[5]))-2)+" l "+2*yscale(pos[5])+" "+3*yscale(pos[5])+" l "+2*yscale(pos[5])+" -"+3*yscale(pos[5])+" z");
				  d3.select("#mask_"+safeNames()[k]+"_"+datatrack+"_"+poskey)
					.select("rect").attr("x", xscale(posx)-yscale(pos[5]));
					break;
				case "pie":
				  var ele = d3.select("#lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+poskey);
				  var s = ele.attr("transform");
				  var t = s.substring(s.indexOf("(")+1, s.indexOf(")"))
									 .split(",");
				  ele.attr("transform", "translate("+xscale(posx)+","+t[1]+")");
				  break;
				case "pie.stack":
				  var ele = d3.select("#lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+poskey);
				  var cls = ele.attr('ref');
				  d3.selectAll("."+cls).attr("transform", function(d){				  	
				  	var y = d3.select(this);
				  	var s= y.attr("transform");
				  	var t = s.substring(s.indexOf("(")+1, s.indexOf(")"))
									 .split(",");
					x.tracklist[trackNames()[k]][datatrack].labpos[y.attr("poskey")] = posx;
					return("translate("+xscale(posx)+"," + t[1] + ")");
				  });
				  
				  break;
				case "dandelion":
				  break;
			 }
			  var txt = d3.select("#lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+poskey);
				  if(txt.size()==1){
				  var s = txt.attr("transform");
				  var t = s.substring(s.indexOf("(")+1, s.indexOf(")"))
									 .split(",");
				  txt.attr("transform", "translate("+xscale(posx)+","+t[1]+")");
			  }
			  d3.select("#lolliplotResizelineL_"+safeNames()[k]+"_"+datatrack+"_"+poskey)
			  	.attr("x1", xscale(posx)-yscale(pos[5]))
			  	.attr("x2", xscale(posx)-yscale(pos[5]));
			  /*d3.select("#lolliplotResizelineR_"+safeNames()[k]+"_"+datatrack+"_"+poskey)
			  	.attr("x1", xscale(posx)+yscale(pos[5]))
			  	.attr("x2", xscale(posx)+yscale(pos[5]));*/
			}
			function draggedResize(){
				var ele = d3.select(this);
				var k = Number(ele.attr("kvalue"));
				var datatrack = ele.attr("datatrack");
				var dx = d3.event.dx;
				var ref = ele.attr("ref");
				if(ref=="+") dx = -1*dx;
				if(dx!=0){
					dx = yscale.invert(yscale(parameter.circleR[safeNames()[k]+"_"+datatrack]) + dx);
					if(dx > 0 && dx < 1){
						parameter.circleR[safeNames()[k]+"_"+datatrack] = dx;
					};
				}
				//renew;
				plotregion.renew();
			}
			function draggedStart(){
				var ele = d3.select(this);
				var k = Number(ele.attr("kvalue"));
				var datatrack = ele.attr("datatrack");
				var tmpstatus = parameter.circleR[safeNames()[k]+"_"+datatrack];
				addNewHistory({
					undo:function(){
						parameter.circleR[safeNames()[k]+"_"+datatrack] = tmpstatus;
						plotregion.renew();
					},
					redo:function(){}
				});
			}
			
			var type = "circle";
			if(typeof(trackdat["stack.factor"])!="undefined"){
				type = "pie.stack";
				// add legend;
				var legend = lolli.append("g");
				var factorLevel = Math.max(...trackdat["stack.factor.order"]);
				for(var i=0; i<factorLevel+1; i++){
					for(var j=0; j<trackdat.start.length; j++){
						if(i==trackdat["stack.factor.order"][j]){
							var col = trackdat.color[Object.keys(trackdat.color)[j]][0];
							var m = textDefaultOptions();
							var id = "legend_"+safeNames()[k]+"_"+datatrack+"_"+col+"_"+trackdat["stack.factor"][j];
							if(typeof(x.markers[id])!="undefined"){
								m = x.markers[id];
								m.vp = legend;
							}else{//init
								m.text = trackdat["stack.factor"][j];
								m.x = xscale.invert(widthF()/2+(i-factorLevel/2)*100);
								m.y = m.fontsize;
								m.anchor = "start";
								m.color = col;
								m.cls = "legend_"+safeNames()[k]+"_"+datatrack;
								m.trackKey = k;
								m.id = id;
								m.angle = 0;
								m.datatrack = datatrack;
								m.ref = m.id;
								m.fill = col;
								m.border = "black";
								m.r = m.fontsize/2;
								m.vp = legend;
							}
							x.markers[m.id] = m;
							var llab = new legendLabel(m);
							break;
						}
					}
				}
			}
			var thisSNP=lolli.append('g');
			
			// add resizable SNPline
			var thisSNPresizeLabel;
			var thisSNPresizeLine = thisSNP.append("line")
							.attr('stroke', 'white')
							.attr('stroke-width', '2px')
							.attr('x1', 0)
							.attr('y1', yscale(pos[3]))
							.attr('x2', widthF())
							.attr('y2', yscale(pos[3]))
							.style("opacity", 0)
							.attr('kvalue', k)
							.attr('datatrack', datatrack)
							.attr('ref', "-")
							.attr("id", "lolliplotResizelineB_"+safeNames()[k]+"_"+datatrack)
							.style("cursor", "ns-resize")
							.on("mouseover", function(){
								d3.select(this).attr("stroke", 'gray').style("opacity", 0.5);
								thisSNPresizeLabel = thisSNP.append("text")
															.attr("fill", "#222222")
															.attr("x", widthF()/2)
															.attr("y", yscale(pos[3])-5)
															.style("font-size", "10px")
															.text("drag to change the height of lollipop stick");
								})
							.on("mouseout", function(){
								d3.select(this).attr("stroke", 'white').style("opacity", 0);
								if(typeof(thisSNPresizeLabel)!="undefined"){
											thisSNPresizeLabel.remove();
										}
								})
							.call(d3.drag().on("drag", function(){
									if(typeof(parameter.YposSTART[safeNames()[k]+"_"+datatrack])=="undefined") parameter.YposSTART[safeNames()[k]+"_"+datatrack] = 0;
									if(d3.event.dy<10) parameter.YposSTART[safeNames()[k]+"_"+datatrack] -= 1-yscale.invert(d3.event.dy/3);
									plotregion.renew();
								}).on("start", function(){
									d3.select(this).style("cursor", "ns-resize").raise().classed("active", true);
									addNewHistory({
										undo:function(){
											parameter.YposSTART[safeNames()[k]+"_"+datatrack] = 0;
											plotregion.renew();
										},
										redo:function(){}
									});
								}).on("end", function(){
									d3.select(this).style("cursor", "default").classed("active", false);
								}));
			
			if(typeof(trackdat.type)!="undefined"){
				if(trackdat.type[0] == "dandelion"){
					thisSNPresizeLine.attr('y1', yscale(pos[1])).attr('y2', yscale(pos[1]));					
					// set group
					var maxgap = (x.end - x.start)/50;
					if(typeof(trackdat.maxgap)!="undefined"){
						maxgap = trackdat.maxgap[0];
					}
					if(typeof(parameter.maxGAP)!="undefined"){
						maxgap *= parameter.maxGAP;
					}
					trackdat.group = [0];
					var groupSize = [1];
					var treeData = [{x:trackdat.start[0],
									 y:groupSize[0],
									 r:1,
									 score:[trackdat.score[0]],
									 mean:trackdat.score[0],
									 data:[1],
									 start:[trackdat.start[0]],
									 color:[trackdat.color[0]||"black"],
									 border:[trackdat.border[0]||"black"],
									 poskey:[0]
									}];
					var method = "mean";
					var colorSets = [];
					if(typeof(trackdat.method)!="undefined"){
						method = trackdat.method[0];
					}
					for(var i=1; i<trackdat.start.length; i++){
						if(typeof(colorSets[trackdat.color[i]])=="undefined"){
							colorSets[trackdat.color[i]]=1;
						}else{
							colorSets[trackdat.color[i]] +=1;
						}
						if(trackdat.start[i] - trackdat.start[i-1]>maxgap){
							trackdat.group[i] = trackdat.group[i-1] + 1;
							groupSize[trackdat.group[i]] = 1;
							treeData[trackdat.group[i]] = {x:0,y:0,r:1,score:[],mean:1,data:[],start:[],color:[],border:[],poskey:[]};
						}else{
							trackdat.group[i] = trackdat.group[i-1];
							groupSize[trackdat.group[i]] = groupSize[trackdat.group[i]] + 1;
							
						}
						treeData[trackdat.group[i]].y = groupSize[trackdat.group[i]];
						treeData[trackdat.group[i]].x = 
							(treeData[trackdat.group[i]].x * treeData[trackdat.group[i]].data.length + trackdat.start[i])/(treeData[trackdat.group[i]].data.length+1);
						treeData[trackdat.group[i]].data.push(1);
						treeData[trackdat.group[i]].r = treeData[trackdat.group[i]].data.length;
						treeData[trackdat.group[i]].score.push(trackdat.score[i]);
						treeData[trackdat.group[i]].start.push(trackdat.start[i]);
						if(typeof(trackdat.border[i])!="undefined"){
							treeData[trackdat.group[i]].border.push(trackdat.border[i]);
						}else{
							treeData[trackdat.group[i]].border.push("black");
						}
						if(typeof(trackdat.color[i])!="undefined"){
							treeData[trackdat.group[i]].color.push(trackdat.color[i]);
						}else{
							treeData[trackdat.group[i]].color.push("black");
						}
						treeData[trackdat.group[i]].poskey.push(i);
					}
					
					var maxScore = 0;
					for(var i=0; i<treeData.length; i++){
						treeData[i].mean = 0;
						switch(method){
							case "count":
								treeData[i].mean = treeData[i].r;
								break;
							case "mean":
								for(var j=0; j<treeData[i].score.length; j++){
									treeData[i].mean += treeData[i].score[j];
								}
								treeData[i].mean = treeData[i].mean/treeData[i].r;
								break;
							default: //mean value of score
								for(var j=0; j<treeData[i].score.length; j++){
									treeData[i].mean += treeData[i].score[j];
								}
								treeData[i].mean = treeData[i].mean/treeData[i].r;
						}
						if(maxScore < treeData[i].mean) maxScore = treeData[i].mean;
					}
						
					var htUnit = (pos[4]-pos[1])/Math.round(...groupSize);
					for(var i=0; i<trackdat.start.length; i++){						
						var snpLine=thisSNP.append('g')
								.attr("class", "lolliplotLine_"+safeNames()[k])
								.attr("kvalue", k)
								.attr("datatrack", datatrack)
								.attr("poskey", i)
								.attr("comp", "lines")
								.on("click", function(){
										var obj=d3.select(this);
										var poskey = Number(obj.attr("poskey"));
										var datatrack = obj.attr("datatrack");
										var k = Number(obj.attr("kvalue"));
										var picked = function(col){
											var tmpstatus = {k:k, 
															 datatrack:datatrack,
															 poskey:poskey,
															 v:clone(x.tracklist[trackNames()[k]][datatrack].border)};
											addNewHistory({
												undo:function(){
													x.tracklist[trackNames()[tmpstatus.k]][tmpstatus.datatrack].border = clone(tmpstatus.v);
													plotregion.renew();
												},
												redo:function(){}
											});
											if(cpCheckAll){
												var thiscol = x.tracklist[trackNames()[k]][datatrack].border[poskey];
												if(typeof(thiscol)=="undefined"){
													if(typeof(x.tracklist[trackNames()[k]][datatrack].border)=="undefined"){
														x.tracklist[trackNames()[k]][datatrack].border = [];
														for(var i=0; i< x.tracklist[trackNames()[k]][datatrack].start.length; i++){
															x.tracklist[trackNames()[k]][datatrack].border[i] = "black";
														}
													}
												}
												for(var j=0; j<x.tracklist[trackNames()[k]][datatrack].border.length; j++){
													if(x.tracklist[trackNames()[k]][datatrack].border[j]==thiscol){
															x.tracklist[trackNames()[k]][datatrack].border[j]=col;
													}
												}
												plotregion.renew();
											}else{
												x.tracklist[trackNames()[k]][datatrack].border[poskey]=col;
												d3.select("#dandelionNodeLinker1_"+safeNames()[k]+"_"+datatrack+"_"+poskey).attr("stroke", col);
											}
										};
										ColorPicker(this, picked, true);
									});
						var bordercolor = "black";
						if(typeof(trackdat.border[i])!="undefined"){
							bordercolor = trackdat.border[i];
						}
						var lastLine = snpLine.append("line")
								 .attr("x1", xscale(trackdat.start[i]))
								 .attr("x2", xscale(trackdat.start[i]))
								 .attr("y1", yscale(pos[0]))
								 .attr("y2", yscale(pos[1]+htUnit*treeData[trackdat.group[i]].mean))
								 .attr("stroke", bordercolor)
								 .attr("id", "dandelionNodeLinker1_"+safeNames()[k]+"_"+datatrack+"_"+i);
					}
					for(var i=0; i<treeData.length; i++){
						treeData[i].mean = 0;
						for(var j=0; j<treeData[i].score.length; j++){
							treeData[i].mean += treeData[i].score[j];
						}
						treeData[i].mean = treeData[i].mean/treeData[i].r;
						
						var nodes = thisSNP.append("g")
										.attr("transform", "translate("+xscale(treeData[i].x)+","+yscale(pos[1]+htUnit*treeData[i].mean)+")");
						var pie = d3.pie().sort(null);
						if(treeData[i].r>10){
							pie.startAngle(Math.PI*.9).endAngle(-Math.PI*.9);
						}else{
							pie.startAngle(Math.PI/2).endAngle(-Math.PI/2);
						}
						nodes.datum(treeData[i].data).selectAll("line")
								.data(pie)
								.enter()
								.append("line")
								.attr("x1", (d,j) => xscale(treeData[i].start[j])-xscale(treeData[i].x))
								.attr("x2", d => Math.sin((d.startAngle+d.endAngle)/2) * -yscale(Math.max(0, 1-(1-pos[5])*Math.sqrt(treeData[i].r)))/2)
								.attr("y1", 0)
								.attr("y2", d => Math.cos((d.startAngle+d.endAngle)/2) * -yscale(Math.max(0, 1-(1-pos[5])*Math.sqrt(treeData[i].r)))/2)
								.attr("stroke", (d,j) => treeData[i].border[j]);
						var circles=nodes.datum(treeData[i].data).selectAll("circle")
								.data(pie)
								.enter()
								.append("circle")
								.attr("r", yscale(pos[5])/4)
								.attr("cx", d => Math.sin((d.startAngle+d.endAngle)/2) * -yscale(Math.max(0, 1-(1-pos[5])*Math.sqrt(treeData[i].r)))/2)
								.attr("cy", d => Math.cos((d.startAngle+d.endAngle)/2) * -yscale(Math.max(0, 1-(1-pos[5])*Math.sqrt(treeData[i].r)))/2)
								.attr("fill", (d,j) => treeData[i].color[j])
								.attr("stroke", (d,j) => treeData[i].border[j])
								.attr("poskey", (d,j) => treeData[i].poskey[j])
								.attr("kvalue", k)
								.attr("datatrack", datatrack)
								.attr("id", (d,j) => "dandelionNode_"+safeNames()[k]+"_"+datatrack+"_"+treeData[i].poskey[j])
								.on("click", function(){
										var obj=d3.select(this);
										var poskey = Number(obj.attr("poskey"));
										var datatrack = obj.attr("datatrack");
										var k = Number(obj.attr("kvalue"));
										var picked = function(col){
											var tmpstatus = {k:k, 
															 datatrack:datatrack,
															 poskey:poskey,
															 v:clone(x.tracklist[trackNames()[k]][datatrack].color)};
											addNewHistory({
												undo:function(){
													x.tracklist[trackNames()[tmpstatus.k]][tmpstatus.datatrack].color = clone(tmpstatus.v);
													plotregion.renew();
												},
												redo:function(){}
											});
											if(cpCheckAll){
												var thiscol = x.tracklist[trackNames()[k]][datatrack].color[poskey];
												if(typeof(thiscol)=="undefined"){
													if(typeof(x.tracklist[trackNames()[k]][datatrack].color)=="undefined"){
														x.tracklist[trackNames()[k]][datatrack].color = [];
														for(var i=0; i< x.tracklist[trackNames()[k]][datatrack].start.length; i++){
															x.tracklist[trackNames()[k]][datatrack].color[i] = "black";
														}
													}
												}
												for(var j=0; j<x.tracklist[trackNames()[k]][datatrack].color.length; j++){
													if(x.tracklist[trackNames()[k]][datatrack].color[j]==thiscol){
															x.tracklist[trackNames()[k]][datatrack].color[j]=col;
													}
												}
												plotregion.renew();
											}else{
												x.tracklist[trackNames()[k]][datatrack].color[poskey]=col;
												d3.select("#dandelionNode_"+safeNames()[k]+"_"+datatrack+"_"+poskey).attr("fill", col);
											}
										};
										cpCheckAll = true;
										ColorPicker(this, picked, true);
									});
						//add resizeLine
					 	nodes.datum(treeData[i].data).selectAll("g")
								.data(pie)
								.enter().append('g')
								.attr("transform", "translate(0,-"+yscale(.5)+")")
								.append("line")
								.attr('stroke', 'white')
								.attr('stroke-width', '2px')
								.attr('x1', d => Math.sin((d.startAngle+d.endAngle)/2) * -yscale(Math.max(0, 1-(1-pos[5])*Math.sqrt(treeData[i].r)))/2 - yscale(pos[5])/4)
								.attr('y1', yscale(0))
								.attr('x2', d => Math.sin((d.startAngle+d.endAngle)/2) * -yscale(Math.max(0, 1-(1-pos[5])*Math.sqrt(treeData[i].r)))/2 - yscale(pos[5])/4)
								.attr('y2', yscale(1))
								.style("opacity", 0)
								.attr('poskey', i)
								.attr('kvalue', k)
								.attr('datatrack', datatrack)
								.attr('ref', "-")
								.attr("id", "lolliplotResizelineL_"+safeNames()[k]+"_"+datatrack+"_"+i)
								.style("cursor", "ew-resize")
								.call(d3.drag().on("drag", draggedResize).on("start", draggedStart));
					}
					//addd legend
					var legend = lolli.append("g");
					for(var i=0; i<Object.keys(colorSets).length; i++){
						var col = Object.keys(colorSets)[i];
						var m = textDefaultOptions();
						if(typeof(x.markers["legend_"+safeNames()[k]+"_"+datatrack+"_"+col])!="undefined"){
							m = x.markers["legend_"+safeNames()[k]+"_"+datatrack+"_"+col];
							m.vp = legend;
						}else{//init
							m.text = "color set " + (i+1);
							m.x = xscale.invert(widthF()/2+(i-Object.keys(colorSets).length/2)*100);
							m.y = m.fontsize;
							m.anchor = "start";
							m.color = col;
							m.cls = "legend_"+safeNames()[k]+"_"+datatrack;
							m.trackKey = k;
							m.id = "legend_"+safeNames()[k]+"_"+datatrack+"_"+col;
							m.angle = 0;
							m.datatrack = datatrack;
							m.ref = m.id;
							m.fill = col;
							m.border = "black";
							m.r = m.fontsize/2;
							m.vp = legend;
						}
						x.markers[m.id] = m;
						var llab = new legendLabel(m);
					}
					return(null);
				}
			}
			
			var setcol1=false;
			var setcol2=false;
			if(typeof(trackdat.color)=="undefined"){
				trackdat.color = [];
				setcol1=true;
			}
			if(typeof(trackdat.color2)=="undefined"){
				trackdat.color2 = [];
				setcol2=true;
			}
			for(var i=0; i<trackdat.start.length; i++){
				var bordercolor = "black";
				if(typeof(trackdat.border[i])!="undefined"){
					bordercolor = trackdat.border[i];
				}
				if(setcol1){
					trackdat.color[i] = bordercolor;
				}
				if(setcol2){
					trackdat.color2[i] = bordercolor;
				}
				var plotSNPline = true;
				if(typeof(trackdat["stack.factor.first"])!="undefined"){
					if(!trackdat["stack.factor.first"][i]){
						plotSNPline = false;
					}
				}
				if(plotSNPline){
					var snpLine=thisSNP.append('g')
							.attr("class", "lolliplotLine_"+safeNames()[k])
							.attr("kvalue", k)
							.attr("datatrack", datatrack)
							.attr("poskey", i)
							.attr("comp", "lines")
							.on("click", function(){
									var obj=d3.select(this);
									var poskey = Number(obj.attr("poskey"));
									var datatrack = obj.attr("datatrack");
									var k = Number(obj.attr("kvalue"));
									var picked = function(col){
										var tmpstatus = {k:k, 
														 datatrack:datatrack,
														 poskey:poskey,
														 v:clone(x.tracklist[trackNames()[k]][datatrack].border)};
										addNewHistory({
											undo:function(){
												x.tracklist[trackNames()[tmpstatus.k]][tmpstatus.datatrack].border = clone(tmpstatus.v);
												plotregion.renew();
											},
											redo:function(){}
										});
										if(cpCheckAll){
											var thiscol = x.tracklist[trackNames()[k]][datatrack].border[poskey];
											if(typeof(thiscol)=="undefined"){
												if(typeof(x.tracklist[trackNames()[k]][datatrack].border)=="undefined"){
													x.tracklist[trackNames()[k]][datatrack].border = [];
													for(var i=0; i< x.tracklist[trackNames()[k]][datatrack].start.length; i++){
														x.tracklist[trackNames()[k]][datatrack].border[i] = "black";
													}
												}
											}
											for(var j=0; j<x.tracklist[trackNames()[k]][datatrack].border.length; j++){
												if(x.tracklist[trackNames()[k]][datatrack].border[j]==thiscol){
														x.tracklist[trackNames()[k]][datatrack].border[j]=col;
												}
											}
											plotregion.renew();
										}else{
											x.tracklist[trackNames()[k]][datatrack].border[poskey]=col;
											var old = d3.select("#lollipopNodeLinkerB_"+safeNames()[k]+"_"+datatrack+"_"+poskey).attr("stroke");
											d3.select("#lollipopNodeLinkerB_"+safeNames()[k]+"_"+datatrack+"_"+poskey).attr("stroke", col);
											d3.select("#lollipopNodeLinker0_"+safeNames()[k]+"_"+datatrack+"_"+poskey).attr("stroke", col);
											d3.select("#lollipopNodeLinker1_"+safeNames()[k]+"_"+datatrack+"_"+poskey).attr("stroke", col);
											d3.select("#lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+poskey).selectAll("circle").attr("stroke", col);
										}
									};
									ColorPicker(this, picked, true);
								});
					snpLine.append("line")
					 .attr("x1", xscale(trackdat.start[i]))
					 .attr("x2", xscale(trackdat.start[i]))
					 .attr("y1", yscale(pos[0]))
					 .attr("y2", yscale(pos[1]))
					 .attr("stroke", bordercolor)
					 .attr("id", "lollipopNodeLinkerB_"+safeNames()[k]+"_"+datatrack+"_"+i);
					snpLine.append("line")
					 .attr("x1", xscale(trackdat.start[i]))
					 .attr("x2", xscale(trackdat.labpos[i]))
					 .attr("y1", yscale(pos[1]))
					 .attr("y2", yscale(pos[2]))
					 .attr("stroke", bordercolor)
					 .attr("id", "lollipopNodeLinker0_"+safeNames()[k]+"_"+datatrack+"_"+i);
					var lastLine = snpLine.append("line")
					 .attr("x1", xscale(trackdat.labpos[i]))
					 .attr("x2", xscale(trackdat.labpos[i]))
					 .attr("y1", yscale(pos[2]))
					 .attr("y2", yscale(pos[4])+yscale(pos[5]))
					 .attr("stroke", bordercolor)
					 .attr("id", "lollipopNodeLinker1_"+safeNames()[k]+"_"+datatrack+"_"+i);
				 
					 //add resizeLine
					 var thisResizeLineL = thisSNP.append('line')
							.attr('stroke', 'white')
							.attr('stroke-width', '2px')
							.attr('x1', xscale(trackdat.labpos[i])-yscale(pos[5]))
							.attr('y1', yscale(0))
							.attr('x2', xscale(trackdat.labpos[i])-yscale(pos[5]))
							.attr('y2', yscale(1))
							.style("opacity", 0)
							.attr('poskey', i)
							.attr('kvalue', k)
							.attr('datatrack', datatrack)
							.attr('ref', "-")
							.attr("id", "lolliplotResizelineL_"+safeNames()[k]+"_"+datatrack+"_"+i)
							.style("cursor", "ew-resize")
							.call(d3.drag().on("drag", draggedResize).on("start", draggedStart));
					 /*thisSNP.append('line')
							.attr('stroke', 'white')
							.attr('stroke-width', '2px')
							.attr('x1', xscale(trackdat.labpos[i])+yscale(pos[5]))
							.attr('y1', yscale(0))
							.attr('x2', xscale(trackdat.labpos[i])+yscale(pos[5]))
							.attr('y2', yscale(1))
							.style("opacity", 0)
							.attr('poskey', i)
							.attr('kvalue', k)
							.attr("ref", "+")
							.attr('datatrack', datatrack)
							.attr("id", "lolliplotResizelineR_"+safeNames()[k]+"_"+datatrack+"_"+i)
							.style("cursor", "ew-resize")
							.call(d3.drag().on("drag", draggedResize));*/
                }
				 //check type
				 if(typeof(trackdat.type)!="undefined"){
				 	type = trackdat.type[i];
				 }
				 
				 switch(type){
				 	case "circle":
				 		{
						var cir=thisSNP.append('g')
							.attr("id", "lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+i)
							.attr("kvalue", k)
							.attr("datatrack", datatrack)
							.attr("poskey", i)
							.attr("comp", "nodes")
							.on("click", function(){
										var obj = d3.select(this);
										var poskey = Number(obj.attr("poskey"));
										var datatrack = obj.attr("datatrack");
										var k = Number(obj.attr("kvalue"));
                                    	var picked = function(col){
                                    		var tmpstatus = {k:k, 
															 datatrack:datatrack,
															 poskey:poskey,
															 v:clone(x.tracklist[trackNames()[k]][datatrack].color)};
											addNewHistory({
												undo:function(){
													x.tracklist[trackNames()[tmpstatus.k]][tmpstatus.datatrack].color = clone(tmpstatus.v);
													plotregion.renew();
												},
												redo:function(){}
											});
											if(cpCheckAll){
												var thiscol = x.tracklist[trackNames()[k]][datatrack].color[poskey];
												if(typeof(thiscol)=="undefined"){
													if(typeof(x.tracklist[trackNames()[k]][datatrack].color)=="undefined"){
														x.tracklist[trackNames()[k]][datatrack].color = [];
														for(var i=0; i< x.tracklist[trackNames()[k]][datatrack].start.length; i++){
															x.tracklist[trackNames()[k]][datatrack].color[i] = "black";
														}
													}
												}
												for(var j=0; j<x.tracklist[trackNames()[k]][datatrack].color.length; j++){
													if(x.tracklist[trackNames()[k]][datatrack].color[j]==thiscol){
															x.tracklist[trackNames()[k]][datatrack].color[j]=col;
													}
												}
												plotregion.renew();
											}else{
												x.tracklist[trackNames()[k]][datatrack].color[poskey]=col;
												plotregion.renew();
                                    		}
                                    	};
                                    	ColorPicker(this, picked, true);
									})
							.call(d3.drag().on("start", dragstarted)
									.on("drag", draggedGroup)
									.on("end", dragended));
						var circenter=yscale(pos[4]);
						var curscore=trackdat.score[i];
						if(curscore==0){
						  cir.append("circle")
								.attr("cx", xscale(trackdat.labpos[i]))
								.attr("cy", circenter)
								.attr("r", yscale(pos[5]))
								.attr("fill", "white")
								.attr("stroke", bordercolor);
						  curscore=1;
						}else{
						  for(var j=0; j<trackdat.score[i]; j++){
							cir.append("circle")
								.attr("cx", xscale(trackdat.labpos[i]))
								.attr("cy", circenter)
								.attr("r", yscale(pos[5]))
								.attr("fill", trackdat.color[i])
								.attr("stroke", bordercolor);
							if(ypos<.5){
								circenter-=2*yscale(pos[5]);
							}else{
								circenter+=2*yscale(pos[5]);
							}
						  }
						}
						//add mask
						var mask = thisSNP.append("defs")
										  .append("mask")
										  .attr("id", "mask_"+safeNames()[k]+"_"+datatrack+"_"+i);
						if(ypos<.5){
							var thisY=yscale(pos[4])+yscale(pos[5])-curscore*2*yscale(pos[5]);
						}else{
							var thisY=yscale(pos[4])-yscale(pos[5]);
						}
						var thisH=curscore*2*yscale(pos[5]);
						mask.append("rect")
							.attr("x", xscale(trackdat.labpos[i])-yscale(pos[5]))
							.attr('y', thisY)
							.attr('width', yscale(pos[5])*2)
							.attr('height', thisH)
							.style("fill", "white")
							.style("opacity", 1);
						/*cir.append("rect")
							.attr("x", xscale(trackdat.labpos[i])-yscale(pos[5]))
							.attr('y', thisY)
							.attr('width', yscale(pos[5])*2)
							.attr('height', thisH)
							.attr("stroke", "black")
							.attr('fill', "white")
							.attr("opacity", .5);*/
						cir.attr("mask", "url(#mask_"+safeNames()[k]+"_"+datatrack+"_"+i+")");
				
						if(typeof(trackdat.textlabel)!="undefined"){
							if(typeof(trackdat.textlabel[i])=="string"){
								var opt = textDefaultOptions();
								opt.id = "lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i;
								opt.angle= -90;
								if(typeof(trackdat["label.parameter.rot"])!="undefined"){
									opt.angle = -trackdat["label.parameter.rot"][i];
								}
								opt.y = circenter;
								opt.x = xscale(trackdat.labpos[i]);
								if(typeof(x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i])!="undefined"){
									opt.x += x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i].x;
									opt.y += x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i].y;
								}
								opt.anchor = "start";
								opt.cls = "nodelabel_"+safeNames()[k];
								opt.trackKey = k;
								opt.text = trackdat.textlabel[i];
								opt.vp = lolli;
								opt.fontsize = x.fontsize["lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i] || opt.fontsize;
								opt.color = x.color["lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i] || opt.color;
								opt.datatrack = datatrack;
								opt.poskey = i;
								opt.colorPickerId = 6;
								var label = new Label(opt);
							}
						}
						}
				 		break;
				 	case "pin":
				 		{
						var circenter=yscale(pos[4]);
						var curscore=trackdat.score[i];
				 		if(ypos<.5){
							var thisY=yscale(pos[4])+yscale(pos[5])-curscore*2*yscale(pos[5]);
						}else{
							var thisY=yscale(pos[4])-yscale(pos[5])+curscore*2*yscale(pos[5]);
						}
				 		lastLine.attr("y2", thisY);
				 		var cir=thisSNP.append('g')
							.attr("id", "lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+i)
							.attr("kvalue", k)
							.attr("datatrack", datatrack)
							.attr("poskey", i)
							.attr("comp", "nodes")
							.call(d3.drag().on("start", dragstarted)
									.on("drag", draggedGroup)
									.on("end", dragended));
				 		
				 		cir.append("circle")
								.attr("cx", xscale(trackdat.labpos[i]))
								.attr("cy", thisY+((ypos<.5?-1:1)*6*yscale(pos[5])))
								.attr("r", 4*yscale(pos[5]))
								.attr("fill", trackdat.color[i])
								.attr("stroke", bordercolor)
								.attr("kvalue", k)
								.attr("datatrack", datatrack)
								.attr("poskey", i)
								.on("click", function(){
											var obj = d3.select(this);
											var poskey = Number(obj.attr("poskey"));
											var datatrack = obj.attr("datatrack");
											var k = Number(obj.attr("kvalue"));
											var picked = function(col){
												var tmpstatus = {k:k, 
																 datatrack:datatrack,
																 poskey:poskey,
																 v:clone(x.tracklist[trackNames()[k]][datatrack].color)};
												addNewHistory({
													undo:function(){
														x.tracklist[trackNames()[tmpstatus.k]][tmpstatus.datatrack].color = clone(tmpstatus.v);
														plotregion.renew();
													},
													redo:function(){}
												});
												if(cpCheckAll){
													var thiscol = x.tracklist[trackNames()[k]][datatrack].color[poskey];
													if(typeof(thiscol)=="undefined"){
														if(typeof(x.tracklist[trackNames()[k]][datatrack].color)=="undefined"){
															x.tracklist[trackNames()[k]][datatrack].color = [];
															for(var i=0; i< x.tracklist[trackNames()[k]][datatrack].start.length; i++){
																x.tracklist[trackNames()[k]][datatrack].color[i] = "black";
															}
														}
													}
													for(var j=0; j<x.tracklist[trackNames()[k]][datatrack].color.length; j++){
														if(x.tracklist[trackNames()[k]][datatrack].color[j]==thiscol){
																x.tracklist[trackNames()[k]][datatrack].color[j]=col;
														}
													}
													plotregion.renew();
												}else{
													x.tracklist[trackNames()[k]][datatrack].color[poskey]=col;
													plotregion.renew();
												}
											};
											ColorPicker(this, picked, true);
										});
						cir.append("circle")
								.attr("cx", xscale(trackdat.labpos[i]))
								.attr("cy", thisY+((ypos<.5?-1:1)*6*yscale(pos[5])))
								.attr("r", 2*yscale(pos[5]))
								.attr("fill", typeof(trackdat.color2[i])=="undefined"?"black":trackdat.color2[i])
								.attr("kvalue", k)
								.attr("datatrack", datatrack)
								.attr("poskey", i)
								.on("click", function(){
											var obj = d3.select(this);
											var poskey = Number(obj.attr("poskey"));
											var datatrack = obj.attr("datatrack");
											var k = Number(obj.attr("kvalue"));
											var picked = function(col){
												var tmpstatus = {k:k, 
																 datatrack:datatrack,
																 poskey:poskey,
																 v:clone(x.tracklist[trackNames()[k]][datatrack].color2)};
												addNewHistory({
													undo:function(){
														x.tracklist[trackNames()[tmpstatus.k]][tmpstatus.datatrack].color2 = clone(tmpstatus.v);
														plotregion.renew();
													},
													redo:function(){}
												});
												if(cpCheckAll){
													var thiscol = x.tracklist[trackNames()[k]][datatrack].color2[poskey];
													if(typeof(thiscol)=="undefined"){
														if(typeof(x.tracklist[trackNames()[k]][datatrack].color2)=="undefined"){
															x.tracklist[trackNames()[k]][datatrack].color2 = [];
															for(var i=0; i< x.tracklist[trackNames()[k]][datatrack].start.length; i++){
																x.tracklist[trackNames()[k]][datatrack].color2[i] = "black";
															}
														}
													}
													for(var j=0; j<x.tracklist[trackNames()[k]][datatrack].color2.length; j++){
														if(x.tracklist[trackNames()[k]][datatrack].color2[j]==thiscol){
																x.tracklist[trackNames()[k]][datatrack].color2[j]=col;
														}
													}
													plotregion.renew();
												}else{
													x.tracklist[trackNames()[k]][datatrack].color2[poskey]=col;
													plotregion.renew();
												}
											};
											ColorPicker(this, picked, true);
										});
				 		cir.append('path')
				 			.attr("fill", trackdat.color[i])
				 			.attr("stroke", "none")
				 			.attr("d", "m "+(xscale(trackdat.labpos[i])-2*yscale(pos[5]))+" "+(thisY+((ypos<.5?-1:1)*2*yscale(pos[5]))-2)+" l "+2*yscale(pos[5])+" "+3*yscale(pos[5])+" l "+2*yscale(pos[5])+" -"+3*yscale(pos[5])+" z");
				 		if(typeof(trackdat.textlabel)!="undefined"){
							if(typeof(trackdat.textlabel[i])=="string"){
								var opt = textDefaultOptions();
								opt.id = "lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i;
								opt.angle= 0;
								if(typeof(trackdat["label.parameter.rot"])!="undefined") opt.angle = -trackdat["label.parameter.rot"][i];
								opt.y = thisY+((ypos<.5?-1:1)*11*yscale(pos[5]));
								opt.x = xscale(trackdat.labpos[i]);
								if(typeof(x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i])!="undefined"){
									opt.x += x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i].x;
									opt.y += x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i].y;
								}
								opt.anchor = "start";
								opt.cls = "nodelabel_"+safeNames()[k];
								opt.trackKey = k;
								opt.text = trackdat.textlabel[i];
								opt.vp = lolli;
								opt.fontsize = x.fontsize["lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i] || opt.fontsize;
								opt.color = x.color["lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i] || bordercolor;
								opt.datatrack = datatrack;
								opt.poskey = i;
								opt.colorPickerId = 6;
								var label = new Label(opt);
							}
						}
						thisResizeLineL.attr("x1", xscale(trackdat.labpos[i])-4*yscale(pos[5]))
						               .attr("x2", xscale(trackdat.labpos[i])-4*yscale(pos[5]));
				 		
				 		}
				 		break;
				 	case "pie":
				 		{
				 		var circenter=yscale(pos[4]);
						var curscore=trackdat.score[i];
				 		if(ypos<.5){
							var thisY=yscale(pos[4])+yscale(pos[5])-curscore*2*yscale(pos[5]);
						}else{
							var thisY=yscale(pos[4])-yscale(pos[5])+curscore*2*yscale(pos[5]);
						}
				 		lastLine.attr("y2", circenter);
				 		var cir=thisSNP.append('g')
							.attr("id", "lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+i)
							.attr("kvalue", k)
							.attr("datatrack", datatrack)
							.attr("poskey", i)
							.attr("comp", "nodes")
							.attr("transform", "translate("+xscale(trackdat.labpos[i])+","+(+circenter+((ypos<.5?-1:1)*6*yscale(pos[5])))+")")
							.on("click", function(){
										var obj = d3.select(this);
										var poskey = Number(obj.attr("poskey"));
										var datatrack = obj.attr("datatrack");
										var k = Number(obj.attr("kvalue"));
                                    	var picked = function(col){
                                    		var tmpstatus = {k:k, 
															 datatrack:datatrack,
															 poskey:poskey,
															 v:clone(x.tracklist[trackNames()[k]][datatrack].color)};
											addNewHistory({
												undo:function(){
													x.tracklist[trackNames()[tmpstatus.k]][tmpstatus.datatrack].color = clone(tmpstatus.v);
													plotregion.renew();
												},
												redo:function(){}
											});
											if(cpCheckAll){
												var thiscol = x.tracklist[trackNames()[k]][datatrack].color[poskey];
												if(typeof(thiscol)=="undefined"){
													if(typeof(x.tracklist[trackNames()[k]][datatrack].color)=="undefined"){
														x.tracklist[trackNames()[k]][datatrack].color = [];
														for(var i=0; i< x.tracklist[trackNames()[k]][datatrack].start.length; i++){
															x.tracklist[trackNames()[k]][datatrack].color[i] = "black";
														}
													}
												}
												for(var j=0; j<x.tracklist[trackNames()[k]][datatrack].color.length; j++){
													if(x.tracklist[trackNames()[k]][datatrack].color[j]==thiscol){
															x.tracklist[trackNames()[k]][datatrack].color[j]=col;
													}
												}
												plotregion.renew();
											}else{
												x.tracklist[trackNames()[k]][datatrack].color[poskey]=col;
												plotregion.renew();
											}
                                    	};
                                    	ColorPicker(this, picked, true);
									})
							.call(d3.drag().on("start", dragstarted)
									.on("drag", draggedGroup)
									.on("end", dragended));
				 		//pie
				 		var data = [trackdat.score[i], Math.max(...trackdat.score)-trackdat.score[i]];
				 		var pie = d3.pie().sort(null);
				 		var arc = d3.arc().innerRadius(0).outerRadius(6*yscale(pos[5]));
				 		function invertColor(hex) {
				 		//https://stackoverflow.com/questions/35969656/how-can-i-generate-the-opposite-color-according-to-current-color
				 			if(typeof(hex)!="string"){
				 				return bordercolor;
				 			}
							if (hex.indexOf('#') === 0) {
								hex = hex.slice(1);
							}
							// convert 3-digit hex to 6-digits.
							if (hex.length === 3) {
								hex = hex[0] + hex[0] + hex[1] + hex[1] + hex[2] + hex[2];
							}
							if (hex.length !== 6) {
								return '#CCCCCC';
							}
							var r = parseInt(hex.slice(0, 2), 16),
								g = parseInt(hex.slice(2, 4), 16),
								b = parseInt(hex.slice(4, 6), 16);
							return (r * 0.299 + g * 0.587 + b * 0.114) > 186
								? '#000000'
								: '#FFFFFF';
						}
				 		var colorSet = [trackdat.color[i], trackdat.color2[i]];
				 		if(colorSet[1]==colorSet[0]){
				 			colorSet[1] = invertColor(colorSet[0]);
				 		}
				 		var getColor = function(i) {
				 			return(colorSet[i]);
				 		}
				 		cir.datum(data).selectAll("path")
				 			.data(pie)
				 			.enter()
				 			.append("path")
				 			.attr("fill", function(d, i) {return getColor(i);} )
				 			.attr("d", arc);
				 		/*cir.append("circle")
								.attr("cx", xscale(trackdat.labpos[i]))
								.attr("cy", circenter+((ypos<.5?-1:1)*6*yscale(pos[5])))
								.attr("r", 4*yscale(pos[5]))
								.attr("fill", trackdat.color[i])
								.attr("stroke", bordercolor);*/
						//label
				 		if(typeof(trackdat.textlabel)!="undefined"){
							if(typeof(trackdat.textlabel[i])=="string"){
								var opt = textDefaultOptions();
								opt.id = "lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i;
								opt.angle= 0;
								if(typeof(trackdat["label.parameter.rot"])!="undefined"){
									opt.angle=-trackdat["label.parameter.rot"][i]
								}
								opt.y = circenter+((ypos<.5?-1:1)*7*yscale(pos[5]));
								opt.x = xscale(trackdat.labpos[i]);
								if(typeof(x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i])!="undefined"){
									opt.x += x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i].x;
									opt.y += x.lolliplotTrackLabel[safeNames()[k]+"_"+datatrack+"_"+i].y;
								}
								opt.anchor = "start";
								opt.cls = "nodelabel_"+safeNames()[k];
								opt.trackKey = k;
								opt.text = trackdat.textlabel[i];
								opt.vp = lolli;
								opt.fontsize = x.fontsize["lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i] || opt.fontsize;
								opt.color = x.color["lolliplotTrackLabel_"+safeNames()[k]+"_"+datatrack+"_"+i] || bordercolor;
								opt.datatrack = datatrack;
								opt.poskey = i;
								opt.colorPickerId = 6;
								var label = new Label(opt);
							}
						}
						thisResizeLineL.attr("x1", xscale(trackdat.labpos[i])-6*yscale(pos[5]))
						               .attr("x2", xscale(trackdat.labpos[i])-6*yscale(pos[5]));
				 		}
				 		break;
				 	case "dandelion":
				 		
				 		break;
				 	case "pie.stack":
				 		var circenter=yscale(pos[4]);
						var curscore=trackdat["stack.factor.order"][i]-1;
				 		if(ypos<.5){
							var thisY=yscale(pos[4])+yscale(pos[5])-curscore*2*yscale(pos[5]);
						}else{
							var thisY=yscale(pos[4])-yscale(pos[5])+curscore*2*yscale(pos[5]);
						}
				 		lastLine.attr("y2", circenter);
				 		var cir=thisSNP.append('g')
							.attr("id", "lolliplotTrackNode_"+safeNames()[k]+"_"+datatrack+"_"+i)
							.attr("kvalue", k)
							.attr("datatrack", datatrack)
							.attr("poskey", i)
							.attr("comp", "nodes")
							.attr("transform", "translate("+xscale(trackdat.labpos[i])+","+(+thisY+((ypos<.5?-1:1)*2*yscale(pos[5])))+")")
							.attr("class", "piestack_"+safeNames()[k]+"_"+datatrack+"_"+trackdat.start[i])
							.attr("ref", "piestack_"+safeNames()[k]+"_"+datatrack+"_"+trackdat.start[i]);
						if(trackdat["stack.factor.first"][i]){
							cir.call(d3.drag().on("start", dragstarted)
									.on("drag", draggedGroup)
									.on("end", dragended));
						}
				 		//pie
				 		var data = [trackdat.score[i], trackdat.score2[i]];
				 		var pie = d3.pie().sort(null);
				 		var arc = d3.arc().innerRadius(0).outerRadius(yscale(pos[5]));
				 		var colorSet = trackdat.color[Object.keys(trackdat.color)[i]];
				 		var getColor = function(i) {
				 			return(colorSet[i]);
				 		}
				 		cir.datum(data).selectAll("path")
				 			.data(pie)
				 			.enter()
				 			.append("path")
				 			.attr("fill", function(d, i) {return getColor(i);} )
				 			.attr("stroke", bordercolor)
				 			.attr("d", arc);
				 		break;
				 }
			}
        };
        
        //gene track
        parameter.geneTrackHeightFactor = [];
        parameter.isoformR=[];
        parameter.trackLayerDataTxt = [];
        var geneTrack = function(layer, track, start, end, xscale, yscale, wscale, line, k, label){
            var self=this;
            self.layer=layer;
            self.redraw = function(){
            	self.layer.selectAll("*").remove();
            	geneTrack(layer, track, start, end, xscale, yscale, wscale, line, k, label);
            };
            var color = track.style.color;
            if(typeof(parameter.geneTrackHeightFactor[trackNames()[k]])=="undefined"){
            	parameter.geneTrackHeightFactor[trackNames()[k]] = 1;
            }
            if(typeof(parameter.isoformR[trackNames()[k]])=="undefined"){
            	parameter.isoformR[trackNames()[k]] = 0;
            }
            if(typeof(color)!="string"){
            	color = color[0];
            }
            var Y=.5;
            if(typeof(track.dat2.start)!="undefined"){
            	if(track.dat2.start.length>0){//plot lolliplot
					self.lolli=layer.append("g").attr("ref", "lollipop");
					lolliplot(self.lolli, track.dat2, xscale, yscale, k, "dat2", .165, "genetrack");
					self.layer=layer.append("g")
							   .attr("group", "genemodel")
							   .attr("transform", "translate(0,"+yscale(.66)+")");
					Y=.25;
            	}
            }
            function draggedResize(){
            	var ele = d3.select(this);
            	if(d3.event.dy> 2){
            		parameter.geneTrackHeightFactor[trackNames()[ele.attr("kvalue")]] *= .75;
            	}
            	if(d3.event.dy< -2){
            		parameter.geneTrackHeightFactor[trackNames()[ele.attr("kvalue")]] *= 1.25;
            	}
            	plotregion.renew();
            }
            function draggedResizeStart(){
            	var k = d3.select(this).attr("kvalue");
            	var tmpstatus = parameter.geneTrackHeightFactor[trackNames()[k]];
            	addNewHistory({
            		undo:function(){
            			parameter.geneTrackHeightFactor[trackNames()[k]] = tmpstatus;
            			plotregion.renew();
            		},
            		redo:function(){}
            	});
            }
            if(track.dat.start.length>0){
                var geneStart=end, geneEnd=start;
                var intronStart=end+5, intronEnd=start-5;
                var fLayer = Math.max(...track.dat.featureLayerID);
                var feature={
                    "utr5" : Y*parameter.geneTrackHeightFactor[trackNames()[k]]/2/fLayer,
                    "CDS"  : Y*parameter.geneTrackHeightFactor[trackNames()[k]]/fLayer,
                    "exon" : Y*parameter.geneTrackHeightFactor[trackNames()[k]]/fLayer,
                    "utr3" : Y*parameter.geneTrackHeightFactor[trackNames()[k]]/2/fLayer,
                    "ncRNA": Y*parameter.geneTrackHeightFactor[trackNames()[k]]/2/fLayer
                };
                for(var j=1; j<=fLayer; j++){
					var data=[];
					var intron=[];
                    
                	var trackLayerData = {start:[],end:[],feature:[],textlabel:label, strand:"*", fill:[], color:color};
                	for(var i=0; i<track.dat.start.length; i++){
                		if(track.dat.featureLayerID[i]==j){
                			trackLayerData.start.push(+track.dat.start[i]);
                			trackLayerData.end.push(+track.dat.end[i]);
                			trackLayerData.feature.push(""+track.dat.feature[i]);
                			if(typeof(track.dat.textlabel)!="undefined") trackLayerData.textlabel = track.dat.textlabel[i] || trackLayerData.textlabel;
                			trackLayerData.strand = track.dat.strand[i];
                			trackLayerData.fill.push(track.dat.fill[i] || color);
                			trackLayerData.color = track.dat.fill[i] || color;
                		}
                	}
					for(var i=0; i<trackLayerData.start.length; i++){
						if(i > 0){
							intronStart = trackLayerData.end[i-1] + 1;
							intronEnd = trackLayerData.start[i] - 1;
							if(intronStart < intronEnd){
								intron.push({
									"x0" : xscale(intronStart) + 5,
									"x1" : xscale(intronEnd)
								});
							}
						}
						var clipx0 = trackLayerData.start[i];
						var clipx1 = trackLayerData.end[i];
						if(clipx0<=end && clipx1>=start){
							if(clipx0<start){
								clipx0 = start;
							}
							if(clipx1>end){
								clipx1 = end;
							}
							data.push({
								"x" : clipx0,
								"h" : feature[trackLayerData.feature[i]],
								"w" : clipx1 - clipx0 + 1,
								"c" : trackLayerData.fill[i]
							});
						}
						if(trackLayerData.start[i]<geneStart){
							geneStart = trackLayerData.start[i];
						}
						if(trackLayerData.end[i]>geneEnd){
							geneEnd = trackLayerData.end[i];
						}
					}
					
					var thisYpos = Y+(j+parameter.isoformR[trackNames()[k]])/(fLayer+1+2*parameter.isoformR[trackNames()[k]])/2;
					var thisY = yscale(thisYpos);
                	// add resize line
                	self.layer.append('line')
				        .attr('stroke', 'white')
				        .attr('stroke-width', '2px')
						.attr("x1", xscale(geneStart))
						.attr("x2", xscale(geneEnd))
						.attr("y1", yscale(thisYpos + Y*parameter.geneTrackHeightFactor[trackNames()[k]]/fLayer/2))
						.attr("y2", yscale(thisYpos + Y*parameter.geneTrackHeightFactor[trackNames()[k]]/fLayer/2))
                        .style("opacity", 0)
                        .attr('fLayer', fLayer)
				        .attr('kvalue', k)
				        .attr("id", "geneTrackResizelineB_"+safeNames()[k]+"_"+fLayer)
				        .style("cursor", "ns-resize")
						.on("mouseover", function(){d3.select(this).attr("stroke", 'gray').style("opacity", 0.5)})
						.on("mouseout", function(){d3.select(this).attr("stroke", 'white').style("opacity", 0)})
                        .call(d3.drag().on("drag", draggedResize).on("start", draggedResizeStart));
					// add a center line
					if(geneStart < start) geneStart = start;
					if(geneEnd > end) geneEnd = end;
					if(geneStart < geneEnd){
						self.layer.append("line")
							 .attr("x1", xscale(geneStart))
							 .attr("x2", xscale(geneEnd))
							 .attr("y1", thisY)
							 .attr("y2", thisY)
							 .attr("stroke", trackLayerData.color)
							 .attr("stroke-width", x.opacity[label]==1 ? "1px" : "10px")
							 .attr("opacity", x.opacity[label])
							 .attr("class", "geneBaseline"+safeNames()[k])
							 .attr("kvalue", k)
							 .attr("ref", j)
							 .on("click", function(){
							 		var obj = d3.select(this);
							 		var featureLayerID = Number(obj.attr("ref"));
							 		var k = Number(obj.attr("kvalue"));
							 		var old = obj.attr("stroke");
							 		var picked = function(col){
							 			var tmpstatus = {k:k, 
														 v:clone(x.tracklist[trackNames()[k]].dat.fill)};
										addNewHistory({
											undo:function(){
												x.tracklist[trackNames()[tmpstatus.k]].dat.fill = clone(tmpstatus.v);
												plotregion.renew();
											},
											redo:function(){}
										});
										for(var i=0; i<x.tracklist[trackNames()[k]].dat.featureLayerID.length; i++){
											if(x.tracklist[trackNames()[k]].dat.featureLayerID[i]==featureLayerID && 
												x.tracklist[trackNames()[k]].dat.fill[i] == old){
												x.tracklist[trackNames()[k]].dat.fill[i] = col;
											}
										}
							 			self.redraw();
							 		};
									ColorPicker(this, picked);
								});
						// add arrows to center line
						plotArrow(self.layer, intron, xscale(start), xscale(end), thisY, trackLayerData.strand, trackLayerData.color, k);
					}
					self.layer.selectAll(".exon"+safeNames()[k])
						.data(data).enter().append("rect")
						.attr("class", "exon track"+safeNames()[k])
						.attr("kvalue", k)
						.attr("x", d=> xscale(d.x))
						.attr("y", d=> yscale(thisYpos + d.h/2))
						.attr("width", d=> wscale(d.w))
						.attr("height", d=> yscale(1 - d.h)>0?yscale(1 - d.h):1)
						.attr("fill", d=> d.c)
						.attr("ref", j)
						.on("click", function(){
									var obj = d3.select(this);
							 		var featureLayerID = Number(obj.attr("ref"));
							 		var k = Number(obj.attr("kvalue"));
							 		var old = obj.attr("fill");
							 		var picked = function(col){
							 			var tmpstatus = {k:k, 
														 v:clone(x.tracklist[trackNames()[k]].dat.fill)};
										addNewHistory({
											undo:function(){
												x.tracklist[trackNames()[tmpstatus.k]].dat.fill = clone(tmpstatus.v);
												plotregion.renew();
											},
											redo:function(){}
										});
										for(var i=0; i<x.tracklist[trackNames()[k]].dat.featureLayerID.length; i++){
											if(x.tracklist[trackNames()[k]].dat.featureLayerID[i]==featureLayerID && 
												x.tracklist[trackNames()[k]].dat.fill[i] == old){
												x.tracklist[trackNames()[k]].dat.fill[i] = col;
											}
										}
							 			self.redraw();
							 		};
									ColorPicker(this, picked);
								});
					// add gene symbols
					var opt = textDefaultOptions();
					opt.angle= typeof(x.rotate[safeNames()[k]+"_"+j])=="undefined"?0:x.rotate[safeNames()[k]+"_"+j];
					opt.y = typeof(x.dataYlabelPos.y[safeNames()[k]+"_0_"+j])=="undefined"?thisY+4:x.dataYlabelPos.y[safeNames()[k]+"_0_"+j];
					opt.x = typeof(x.dataYlabelPos.x[safeNames()[k]+"_0_"+j])=="undefined"?xscale(geneStart) - 10:xscale(x.dataYlabelPos.x[safeNames()[k]+"_0_"+j]);
					opt.anchor = "end";
					opt.cls = "trackLayerLabel_"+safeNames()[k];
					opt.trackKey = k;
					opt.text = trackLayerData.textlabel;
					if(typeof(parameter.trackLayerDataTxt[safeNames()[k]+"_"+j])!="undefined") opt.text = parameter.trackLayerDataTxt[safeNames()[k]+"_"+j];
					parameter.trackLayerDataTxt[safeNames()[k]+"_"+j] = opt.text;
					opt.vp = self.layer;
					opt.fontsize = x.fontsize[safeNames()[k]+"_"+j] || defaultFontSize;
					opt.color = x.color[trackNames()[k]] || opt.color;
					opt.colorPickerId = 3;
					opt.ref = j;
					self.ylabel = new Label(opt);
				}
            }
            return(self);
        };
        
        //lollipop plot track
        var lolliTrack = function(layer, track, start, end, xscale, yscale, wscale, line, k, label){
        	var self = this;
            var color = track.style.color;
            if(typeof(color)=="string"){
            	color = [color, color];
            }
            var ypos = [0.45, 0.55];
            if(typeof(track.dat.start)=="undefined"){
                ypos[1]=.825;
            }
            if(typeof(track.dat2.start)=="undefined"){
                ypos[0]=.165;
            }
            
            if(typeof(track.dat.start)!="undefined"){
            	if(track.dat.start.length>0){//plot lolliplot
					self.lolli1=layer.append("g").attr("group", "lollipop");
					// add a center line
					self.lolli1.append("line")
						 .attr("x1", xscale(start))
						 .attr("x2", xscale(end))
						 .attr("y1", yscale(ypos[0]))
						 .attr("y2", yscale(ypos[0]))
						 .attr("stroke", color[0])
						 .attr("stroke-width", x.opacity[label]==1 ? "1px" : "10px")
						 .attr("opacity", x.opacity[label])
						 .attr("id", "geneBaseline_"+safeNames()[k]+"_dat");
					self.lolli1.append("line")
						 .attr("x1", xscale(start))
						 .attr("x2", xscale(end))
						 .attr("y1", yscale(ypos[0]))
						 .attr("y2", yscale(ypos[0]))
						 .attr("stroke", color[0])
						 .attr("stroke-width", "10px")
						 .attr("opacity", 0)
						 .attr("class", "geneBaseline_"+safeNames()[k])
						 .attr("kvalue", k)
						 .on("click", function(){
						 	var target = d3.select(this);
						 	var picked = function(col){
								var k = Number(target.attr("kvalue"));
								var tmpstatus = {k:k,
												v:x.tracklist[trackNames()[k]].style.color[0]};
								addNewHistory({
									undo:function(){
										x.tracklist[trackNames()[tmpstatus.k]].style.color[0] = tmpstatus.v;
										plotregion.renew();
									},
									redo:function(){}
								});
								x.tracklist[trackNames()[k]].style.color[0]=col;
								d3.select("#geneBaseline_"+safeNames()[k]+"_dat").attr("stroke", col);
						 	};
						 	ColorPicker(this, picked);
						 });
					self.lolliplot1=lolliplot(self.lolli1, track.dat, xscale, yscale, k, "dat", ypos[0], "lollipopData");
            	}
            }
            
            if(typeof(track.dat2.start)!="undefined"){
            	if(track.dat2.start.length>0){//plot lolliplot
					self.lolli2=layer.append("g").attr("group", "lollipop");
					// add a center line
					self.lolli2.append("line")
						 .attr("x1", xscale(start))
						 .attr("x2", xscale(end))
						 .attr("y1", yscale(ypos[1]))
						 .attr("y2", yscale(ypos[1]))
						 .attr("stroke", color[1])
						 .attr("stroke-width", x.opacity[label]==1 ? "1px" : "10px")
						 .attr("opacity", x.opacity[label])
						 .attr("id", "geneBaseline_"+safeNames()[k]+"_dat2");
					self.lolli2.append("line")
						 .attr("x1", xscale(start))
						 .attr("x2", xscale(end))
						 .attr("y1", yscale(ypos[1]))
						 .attr("y2", yscale(ypos[1]))
						 .attr("stroke", color[1])
						 .attr("stroke-width", "10px")
						 .attr("opacity", 0)
						 .attr("class", "geneBaseline_"+safeNames()[k])
						 .attr("kvalue", k)
						 .on("click", function(){
						 	var target = d3.select(this);
						 	var picked = function(col){
								var k = Number(target.attr("kvalue"));
								var tmpstatus = {k:k,
												v:x.tracklist[trackNames()[k]].style.color[1]};
								addNewHistory({
									undo:function(){
										x.tracklist[trackNames()[tmpstatus.k]].style.color[1] = tmpstatus.v;
										plotregion.renew();
									},
									redo:function(){}
								});
								x.tracklist[trackNames()[k]].style.color[1]=col;
								d3.select("#geneBaseline_"+safeNames()[k]+"_dat2").attr("stroke", col);
						 	};
						 	ColorPicker(this, picked);
						 });
					self.lolliplot2=lolliplot(self.lolli2, track.dat2, xscale, yscale, k, "dat2", ypos[1], "lollipopData");
            	}
            }
            // add name symbols
            var opt = textDefaultOptions();
            opt.angle= typeof(x.rotate[trackNames()[k]])=="undefined"?0:x.rotate[trackNames()[k]];
			opt.y = typeof(x.dataYlabelPos.y[safeNames()[k]+"_0_text"])=="undefined"?yscale(0.5)+4:x.dataYlabelPos.y[safeNames()[k]+"_0_text"];
			opt.x = typeof(x.dataYlabelPos.x[safeNames()[k]+"_0_text"])=="undefined"?xscale(start) - 10:xscale(x.dataYlabelPos.x[safeNames()[k]+"_0_text"]);
            opt.anchor = "end";
            opt.cls = "dataYlabel_"+safeNames()[k];
            opt.trackKey = k;
            opt.text = label;
            opt.vp = layer;
            opt.fontsize = x.fontsize[label];
            opt.color = x.color[trackNames()[k]] || color[0];
            opt.colorPickerId = 3;
            self.ylabel = new Label(opt);
            
            return(self);
        };
        
        //data track
        var dataTrack = function(layer, track, start, end, xscale, yscale, line, k, height){
        	var self = this;
            var color = track.style.color;
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
                self.dat = layer.append("path")
                    .datum(data)
                    .attr("fill", color[0])
                    .attr("stroke", "none")
                    .attr("d", line)
                    .attr("class", "dataPath1 track"+safeNames()[k])
                    .attr("kvalue", k)
                    .on("click", function(){
                    	var obj = d3.select(this);
                    	var k=obj.attr("kvalue");
                    	var picked = function(col){
                    		var tmpstatus = {k:k, 
											 v:x.tracklist[trackNames()[k]].style.color[0]};
							addNewHistory({
								undo:function(){
									x.tracklist[trackNames()[tmpstatus.k]].style.color[0] = tmpstatus.v;
									plotregion.renew();
								},
								redo:function(){}
							});
                    		x.tracklist[trackNames()[k]].style.color[0] = col;
                    		obj.attr("fill", col);
                    		d3.select(".dataBaseline"+safeNames()[k]).attr("stroke", col);
                    	};
                    	ColorPicker(this, picked);
                    });
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
                self.dat2 = layer.append("path")
                    .datum(data)
                    .attr("fill", color[1])
                    .attr("stroke", "none")
                    .attr("d", line)
                    .attr("class", "dataPath2 track"+safeNames()[k])
                    .attr("kvalue", k)
                    .on("click", function(){
                    	var obj = d3.select(this);
                    	var k=obj.attr("kvalue");
                    	var picked = function(col){
                    		var tmpstatus = {k:k, 
											 v:x.tracklist[trackNames()[k]].style.color[1]};
							addNewHistory({
								undo:function(){
									x.tracklist[trackNames()[tmpstatus.k]].style.color[1] = tmpstatus.v;
									plotregion.renew();
								},
								redo:function(){}
							});
                    		x.tracklist[trackNames()[k]].style.color[1] = col;
                    		obj.attr("fill", col);
                    	};
                    	ColorPicker(this, picked);
                    });
            }
            // add a line at 0
            self.baseline = layer.append("line")
                 .attr("x1", xscale(start))
                 .attr("x2", xscale(end))
                 .attr("y1", yscale(0))
                 .attr("y2", yscale(0))
                 .attr("stroke", color[0])
                 .attr("stroke-width", x.opacity[x.name[k]]==1 ? "1px" : "10px")
                 .attr("opacity", x.opacity[x.name[k]])
                 .attr("class", "dataBaseline"+safeNames()[k])
                 .attr("kvalue", k);
            
            self.yaxis = new yaxis(layer, yscale, track.ylim, height, trackNames()[k], k);
            
            return(self);
        };
        
        
        // xaxis
        parameter.xaxOpt = {angle: 0, fontsize: defaultFontSize, ddx: 0, ddy: 0, color: "#000"};
        var xaxis = function(layer){
        	var self=this;
        	self.ticks=[];
            self.g=layer.append("g").attr("class", "xaxis")
						.attr("transform", "translate(0," + heightF() + ")")
						.call(d3.axisBottom(xscale()).tickSize(5).ticks(5));
            d3.selectAll(".xaxis>.tick>text").each(function(d, i){
            	var thisTick=d3.select(this);
            	var vp = d3.select(this.parentNode);
            	var opt = textDefaultOptions();
            	opt.id = "xaxis_"+i;
            	opt.text = thisTick.text();
            	opt.color = parameter.xaxOpt.color;
            	opt.x = Number(thisTick.attr("x"))+parameter.xaxOpt.ddx;
            	opt.y = Number(thisTick.attr("y"))+parameter.xaxOpt.ddy;
            	opt.dy = thisTick.attr("dy");
            	opt.cls = "xaxis_tick";
            	opt.angle = parameter.xaxOpt.angle;
            	opt.fontsize = parameter.xaxOpt.fontsize;
            	opt.vp = vp;
            	self.ticks.push(new Label(opt));
            	thisTick.remove();
            });
            self.remove = function(){
            	self.g.remove();
            }
            return(self);
        };
        
        // yaxis
        var yaxis = function(layer, yscale, ylim, height, label, k){
        	var self = this;
        	self.g = layer.append("g");
        	self.ticks = [];
        	// yaxis
        	var yaxisStyle = x.tracklist[trackNames()[k]].style;
        	if(yaxisStyle.yaxis.draw){
        		if(yaxisStyle.yaxis.main){
					self.tickG = self.g.append("g").attr("class", "yaxisG_"+safeNames()[k])
						.call(d3.axisLeft(yscale).tickSize(5).tickValues(ylim));
				}else{
					self.tickG = self.g.append("g").attr("class", "yaxisG_"+safeNames()[k]).attr("transform", "translate("+widthF()+",0)")
						.call(d3.axisRight(yscale).tickSize(5).tickValues(ylim));
				}
				d3.selectAll(".yaxisG_"+safeNames()[k]+">.tick>text").each(function(d, i){
					var thisTick=d3.select(this);
					var vp = d3.select(this.parentNode);
					var opt = textDefaultOptions();
					opt.id = "yaxis_"+safeNames()[k]+"_"+i;
					opt.text = thisTick.text();
					opt.y = typeof(x.dataYlabelPos.y[safeNames()[k]+"_"+i+"_text_tick"])=="undefined"?thisTick.attr("y"):x.dataYlabelPos.y[safeNames()[k]+"_"+i+"_text_tick"];
					opt.x = typeof(x.dataYlabelPos.x[safeNames()[k]+"_"+i+"_text_tick"])=="undefined"?thisTick.attr("x"):x.dataYlabelPos.x[safeNames()[k]+"_"+i+"_text_tick"];
					opt.dy = thisTick.attr("dy");
					opt.anchor = yaxisStyle.yaxis.main?"end":"start";
					opt.cls = "yaxis_tick_"+safeNames()[k];
					opt.trackKey = k;
					opt.datatrack = i;
					opt.fontsize = x.fontsize[safeNames()[k]+"__tick"]||defaultFontSize;
					opt.vp = vp;
					self.ticks.push(new Label(opt));
					thisTick.remove();
				});
            }
            self.remove = function(){
            	self.g.remove();
            };
            // trackName
            var opt = textDefaultOptions();
            opt.angle= typeof(x.rotate[trackNames()[k]])=="undefined"?(yaxisStyle.yaxis.main?-90:0):x.rotate[trackNames()[k]];
			opt.y = typeof(x.dataYlabelPos.y[safeNames()[k]+"_0_text"])=="undefined"?height/2:x.dataYlabelPos.y[safeNames()[k]+"_0_text"];
			opt.x = typeof(x.dataYlabelPos.x[safeNames()[k]+"_0_text"])=="undefined"?-8:xscale()(x.dataYlabelPos.x[safeNames()[k]+"_0_text"]);
            opt.anchor = yaxisStyle.yaxis.main?"middle":"end";
            opt.cls = "dataYlabel_"+safeNames()[k];
            opt.trackKey = k;
            opt.text = label;
            opt.vp = self.g;
            opt.fontsize = x.fontsize[label];
            opt.color = x.color[trackNames()[k]] || yaxisStyle.ylabgp.col || opt.color;
            opt.colorPickerId = 3;
            self.ylabel = new Label(opt);
            return(self);
        };
        
        // get real time width and height
        var widthF = function() {
            return(+svg.attr("width") - margin.left - margin.right);
        };
        var heightF = function() {
            return(+svg.attr("height") - margin.top - margin.bottom);
        };
        
        var Track = function(plotregion, k, currH, vspace){
        	var self=this;
        	self.bottom=currH;
        	self.height = function(){
        		return(x.height[trackNames()[k]]*heightF() - 2*vspace);
        	};
        	self.track = plotregion.g.append("g")
						  .attr("kvalue", k)
						  .attr("id", "track"+safeNames()[k])
						  .attr("type", x.type[trackNames()[k]])
						  .attr("width", widthF())
						  .attr("height", self.height())
						  .attr("transform", "translate(0," + (currH + vspace) +")");
			self.remove = function(){
				self.track.remove();
			};
            //make track drag-able
            //move tracks order
            var currentId;
			function dragstarted(d) {
			  d3.select(this).style("cursor", "move").raise().classed("active", true);
			  var tmpstatus = clone(x.name);
			  addNewHistory({
								undo:function(){
									x.name=tmpstatus;
									plotregion.renew();
								},
								redo:function(){}
							  });
			  currentId = Number(d3.select(this).attr("kvalue"));
			}
			function dragended(d) {
			  d3.select(this).style("cursor", "default").classed("active", false);
			}
			//change track height
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
				plotregion.renew();
			};
			function draggedTrack(d){
				var tk = d3.select(this);
				tk.attr("transform", "translate(0," + d3.event.y + ")");
				checkTrackHeight(d3.event.y);
			}
            self.track.call(d3.drag().on("start", dragstarted)
                                    .on("drag", draggedTrack)
                                    .on("end", dragended));
            //adjustable track height
            self.coor = [0,0];
            self.trackBottom = self.track.append("line")
								   .attr("kvalue", k)
								   .attr("id", "trackbottom"+safeNames()[k])
								   .attr("x1", 0)
								   .attr("y1", self.height())
								   .attr("x2", widthF())
								   .attr("y2", self.height())
								   .attr("stroke", "white")
								   .attr("stroke-width", "3px")
								   .attr("opacity", 0)
								   .style("cursor", "ns-resize")
									.on("mouseover", function(){
										d3.select(this).attr("stroke", 'gray').style("opacity", 0.5);
										self.trackBottomText = self.track.append("text")
																	.attr("fill", "#222222")
																	.attr("x", widthF()/2)
																	.attr("y", self.height()-5)
																	.style("font-size", "10px")
																	.text("drag to change the height of track");
									})
									.on("mouseout", function(){
										d3.select(this).attr("stroke", 'white').style("opacity", 0);
										if(typeof(self.trackBottomText)!="undefined"){
											self.trackBottomText.remove();
										}
									})
								   .call(d3.drag()
								   .on("start", function(d){
									   self.coor = d3.mouse(svg.node());
									   var tmpstatus = clone(x.height);
									   addNewHistory({
										undo:function(){
											x.height=tmpstatus;
											plotregion.renew();
										},
										redo:function(){}
									  });
								   })
								   .on("drag", function(d){
									   var obj = d3.select(this);
									   var dy = d3.mouse(svg.node())[1] - self.coor[1];
									   self.coor = d3.mouse(svg.node());
									   var k = +obj.attr("kvalue");
									   var thisH = +obj.attr("y1") + dy;
									   if(thisH > 0){
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
										   plotregion.renew();
									   }
								   }));
			var currTrack = x.tracklist[trackNames()[k]];
			self.yscale = function(){
				return(d3.scaleLinear().domain(currTrack.ylim).rangeRound([self.height(), 0]));
			}
			//line x y;
			var line = d3.line()
				.x(function(d) {return xscale()(d.x);})
				.y(function(d) {return self.yscale()(d.y);});
			switch(x.type[trackNames()[k]]){
				case "data":
					self.plot = new dataTrack(self.track, currTrack, x.start, x.end, xscale(), self.yscale(), line, k, self.height());
					break;
				case "lollipopData":
					lolliTrack(self.track, currTrack, x.start, x.end, xscale(), self.yscale(), wscale(), line, k, trackNames()[k]);
					break;
				case "gene":
				case "transcript":
					self.plot = new geneTrack(self.track, currTrack, x.start, x.end, xscale(), self.yscale(), wscale(), line, k, trackNames()[k]);
					break;
				default:
					console.log(x.type[trackNames()[k]]);
			}
				
        	return(self);
        }
        
        var wscale = function(){
			return(d3.scaleLinear().domain([0, x.end-x.start+1]).rangeRound([0, widthF()]));
        };
		var xscale = function(){
			return(d3.scaleLinear().domain([x.start, x.end]).rangeRound([0, widthF()]));
        };
            
        var Draw = function(vspace=10){
        	var self=this;
            // svg > g
            self.g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
            
            //svg > g > rect: plot region. click to remove highlight.
            /*self.plotRegion=self.g.append("rect")
				 .attr("x", -margin.left)
				 .attr("y", margin.top)
				 .attr("width", svg.attr("width"))
				 .attr("height", svg.attr("height"))
				 .attr("fill", "white")
				 .attr("stroke", "none")
				 .style("opacity", 0)
				 .on("click", function(){
				 //click to remove highligth
				});*/
            
			self.tracks = [];
            self.plotTracks = function(){
            	var currH = 0;
            	self.tracks = [];
				//plot track
				for(var k=0; k<trackNames().length; k++){
					self.tracks.push(new Track(self, k, currH, vspace));
					currH +=x.height[trackNames()[k]]*heightF();
				}
            };
            self.plotTracks();
            self.renewTracks = function(){
            	for(var i=0; i<self.tracks.length; i++){
            		self.tracks[i].remove();
            	}
            	self.plotTracks();            	
            };
            
            //xaxis
            self.x_axis = new xaxis(self.g);
            self.renewXaxis = function(){
            	self.x_axis.remove();
            	self.x_axis = new xaxis(self.g);
            };
            // resize svg size
            var resizeBtnSize = 20;
            self.resizeBtn = svg.append("rect").attr("fill", "white")
                               .attr("x", svg.attr("width")-resizeBtnSize)
                               .attr("y", svg.attr("height")-resizeBtnSize)
                               .attr("width", resizeBtnSize)
                               .attr("height", resizeBtnSize)
                               .attr("opacity", 0)
                               .style("cursor", "nwse-resize")
                               .call(d3.drag().on("drag", function(d){
									svg.attr("width", d3.event.x).attr("height", d3.event.y);
									/*self.plotRegion.attr("width", svg.attr("width"))
												 .attr("height", svg.attr("height"));*/
									self.renew();
									if(typeof(self.ruler)!="undefined") ruler.rulemove(d3.mouse(this));
                               }).on("start", function(d){
									   var tmpstatus = {w:svg.attr("width"),h:svg.attr("height")};
									   addNewHistory({
										undo:function(){
											svg.attr("width", tmpstatus.w).attr("height", tmpstatus.h);
											self.renew();
										},
										redo:function(){}
									  });
								   }));
            self.renewResizeBtn = function(){
            	self.resizeBtn.attr("x", svg.attr("width")-resizeBtnSize)
							.attr("y", svg.attr("height")-resizeBtnSize);
            };
            
            self.marginBox = new Margin(self);
            self.renewMarginBox = function(){
				if(typeof(self.marginBox)!="undefined"){
					self.marginBox.remove();
				}
				self.marginBox = new Margin(self);
            };
			self.renew = function(){
				self.renewXaxis();
				self.renewMarginBox();
				self.renewResizeBtn();
				self.renewTracks();
				self.markers.redraw();
			};					
            self.markers = new Marker();
            
            //context menu
            self.menu = contextMenu(self.markers);
            svg.on('contextmenu', function(){
            	d3.event.preventDefault();
            	self.menu(d3.mouse(this)[0], d3.mouse(this)[1]);
            });
            return(self);
        };
        
        //clear canvas
        svg.selectAll("*").remove();
        var plotregion = Draw();
        var ruler = new Ruler();
        console.log(x);
        
        // set undo
        function clone(obj) {//from https://stackoverflow.com/questions/728360/how-do-i-correctly-clone-a-javascript-object
			var copy;
			// Handle the 3 simple types, and null or undefined
			if (null == obj || "object" != typeof obj) return obj;

			// Handle Date
			if (obj instanceof Date) {
				copy = new Date();
				copy.setTime(obj.getTime());
				return copy;
			}

			// Handle Array
			if (obj instanceof Array) {
				copy = [];
				for (var i = 0, len = obj.length; i < len; i++) {
					copy[i] = clone(obj[i]);
				}
				return copy;
			}

			// Handle Object
			if (obj instanceof Object) {
				copy = {};
				for (var attr in obj) {
					if (obj.hasOwnProperty(attr)) copy[attr] = clone(obj[attr]);
				}
				return copy;
			}

			throw new Error("Unable to copy obj! Its type isn't supported.");
		}
        var history = [];
        var historyIndex = 0;
        function addNewHistory(item){
        	if(history.length>0) history=history.slice(0, historyIndex);
        	history.push(item);
			historyIndex++;
        }
        function undo(){
        	if(history.length>0 && historyIndex > 0){
        		historyIndex--;
        		var action = history[historyIndex];
        		if(typeof(action)!="undefined") action.undo();
        	}
        }
        function redo(){
        	if(history.length>historyIndex){
        		historyIndex++;
        		var action = history[historyIndex];
        		action.redo();
        	}
        }
        d3.select("#undo").on("click", undo);
      },

      resize: function(width, height) {
        svg.attr("width", width)
           .attr("height", height);
      },

      svg: svg
    };
  }
});
