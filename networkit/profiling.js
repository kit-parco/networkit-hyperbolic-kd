/*
	file: profiling.js
	author: Mark Erb
		
	This file contains all JavaScript scripts for the profiling layout 
	
	The data will be automatically embedded as singelton into the original
	IPython Layout 
	
		html > head > script id="NetworKit_script"
		
	In addition a function for hiding the Overlay is defined in 
	"profiling.py"
	
	To prevent conflicts with the IPython Notebook, prove that every function
	begins with	"NetworKit_".
*/

/*
	Markup the data from "profiling.profile.html"
	
	Arguments:
		id: ID of the "NetworKit_Page" to markup
 */
function NetworKit_pageEmbed(id)
{
	var elements = document.getElementById(id).getElementsByClassName("Plot");
	var i, j;
	for (i=0; i<elements.length; i++) {
		elements[i].id = id + "_Plot_" + i;
		var data = elements[i].getAttribute("data-image").split("|");
		elements[i].removeAttribute("data-image");
		var content = 
			"<div class=\\"Image\\" id=\\"" + elements[i].id + "_Image\\">" +
			"  <div class=\\"Title\\">" + elements[i].title + "</div>" +
			"</div>";
		elements[i].innerHTML = content;
		elements[i].setAttribute("data-image-index", 0);
		elements[i].setAttribute("data-image-length", data.length);
		for (j=0; j<data.length; j++) {
			elements[i].setAttribute("data-image-" + j, data[j]);
		}
		NetworKit_plotUpdate(elements[i]);
		elements[i].onclick = function (e) {
			NetworKit_overlayShow((e.target) ? e.target : e.srcElement);
		}
	}
}


/*
	Update Thumbnails.
	
	Arguments:
		source: Element (class=NetworKit_Plot) to update
*/
function NetworKit_plotUpdate(source)
{
	var index = source.getAttribute("data-image-index");
	var data = source.getAttribute("data-image-" + index);
	var image = document.getElementById(source.id + "_Image");
	image.style.backgroundImage = "url(" + data + ")";
}


/*
	"Hide/show" an element.
	
	Arguments:
		id: ID of the element to hide/show
		show: true  -> show
			  false -> hide
*/
function NetworKit_showElement(id, show)
{
	var element = document.getElementById(id);
	element.style.display = (show) ? "block" : "none";
}


/*
	Show Overlay.
	
	Arguments:
		source: element to get data for the overlay
*/
function NetworKit_overlayShow(source)
{
	NetworKit_overlayUpdate(source);
	NetworKit_showElement("NetworKit_Overlay", true);
}


/*
	Update data of overlay.

	Arguments:
		source: element to get data for the overlay
*/
function NetworKit_overlayUpdate(source)
{
	document.getElementById("NetworKit_Overlay_Title").innerHTML = source.title;
	var index = source.getAttribute("data-image-index");
	var data = source.getAttribute("data-image-" + index);
	var image = document.getElementById("NetworKit_Overlay_Image");
	image.setAttribute("data-id", source.id);
	image.style.backgroundImage = "url(" + data + ")";
	var link = document.getElementById("NetworKit_Overlay_Toolbar_Bottom_Save");
	link.href = data;
	link.download = source.title + ".svg";
}


/*
	Shift through the possible images of the overlay.
	
	Arguments:
		delta: >0 next delta-th picture
		       <0 previous delta-th picture
*/
function NetworKit_overlayImageShift(delta)
{
	var image = document.getElementById("NetworKit_Overlay_Image");
	var source = document.getElementById(image.getAttribute("data-id"));
	var index = parseInt(source.getAttribute("data-image-index"));
	var length = parseInt(source.getAttribute("data-image-length"));
	var index = (index+delta) % length;
	if (index < 0) {
		index = length + index;
	}
	source.setAttribute("data-image-index", index);
	NetworKit_overlayUpdate(source);
}