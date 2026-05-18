
var cell = 0;
var channels = 3;

macro "Cell Selector Tool -C44f-o4499" {
	cell += 1;
	showStatus("Changed cell number to " + cell);
	
	getCursorLoc(x, y, t, flags);
	addResult("Start", cell, x, y, t/channels+1);
}

macro "Macro 1 [a]" { 
	getCursorLoc(x, y, t, flags);
	addResult("Bud", cell, x, y, t/channels+1);
} 

macro "Macro 1 [c]" { 
	getCursorLoc(x, y, t, flags);
	addResult("Cytokinesis", cell, x, y, t/channels+1);
} 

macro "Macro 1 [d]" { 
	getCursorLoc(x, y, t, flags);
	addResult("Death", cell, x, y, t/channels+1);

} 

macro "Macro 1 [w]" {
	getCursorLoc(x, y, t, flags);
	addResult("Washed", cell, x, y, t/channels+1);
} 


macro "Macro 1 [n]" {
	cell += 1;
	showStatus("Changed cell number to " + cell);

}

macro "Macro 1 [p]" {
	cell -= 1;
	showStatus("Changed cell number to " + cell);
}

macro "Macro 1 [s]" {
	Dialog.create("Cell Annotator");
	Dialog.addNumber("Number of channels", channels);
	Dialog.addNumber("Current cell number", cell);
	Dialog.show();
	
	channels = Dialog.getNumber();
	cell = Dialog.getNumber();
}

macro "Macro 1 [e]" {
	// export table to correct format
	for (row = 0; row < nResults; row++) {
		event = getResultString("Event", row);
		cell = getResult("Cell", row);
		x = getResult("x", row);
		y = getResult("y", row);
		t = getResult("t", row);

		if (event == "Start") {
			print("Cell" + cell, "Start:", t, "Cell position:", x, y);
		}
		else if (event == "Bud") {
			print("Bud:", t);
		}
		else if (event == "Cytokinesis") {
			print("Cytokinesis:", t);
		}
		else if (event == "Death") {
			print("Death:", t);
		}
		else if (event == "Washed") {
			print("Washed:", t);
		}
	}
}


function addResult(event, cell, x, y, t) {
	row = nResults;
	setResult("Event", row, event);
	setResult("Cell", row, cell);
	setResult("x", row, x);
	setResult("y", row, y);
	setResult("t", row, t);

	Table.sort("t");
	Table.sort("Cell");
	
	updateResults;
}


