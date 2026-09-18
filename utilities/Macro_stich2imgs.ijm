gap = 20;

titles = getList("image.titles");

selectWindow(titles[0]);
w1 = getWidth();
h1 = getHeight();

selectWindow(titles[1]);
w2 = getWidth();
h2 = getHeight();
run("Select All");
run("Copy");

newImage("stitched", "8-bit black", w1 + gap + w2, h1, 1);

makeRectangle(0, 0, w2, h2);
run("Paste");

selectWindow(titles[0]);
run("Select All");
run("Copy");

selectWindow("stitched");
makeRectangle(w2 + gap, 0, w1, h1);
run("Paste");

run("Select None");