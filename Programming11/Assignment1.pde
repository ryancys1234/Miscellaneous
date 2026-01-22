size(780,500);
background(178,185,255);
noStroke();
//ground
fill(85,196,93);
rect(0,400,780,100);
//sun
fill(146,255,5);
arc(600,180,125,125, radians(100), radians(420));
//mountains
fill(38,42,209);
triangle(0,400,200,100,450,400);
triangle(300,400,600,180,780,400);
//robot
fill(46,134,126);
arc(250,430,50,50,0,PI);

stroke(46,134,126);
strokeWeight(2);

line(310,90,275,197);
line(250,145,275,197);
line(350,210,275,197);

line(260,210,260,350);
line(273,210,273,350);
line(286,210,286,350);

noStroke();
fill(0,0,0);
rect(200,310,100,120);

ellipse(275,197,98,98);
fill(#FFFFFF);
ellipse(285,190,30,30);
fill(#82D1D3);
ellipse(298,172,7.3,7.3);
ellipse(255,187,9,9);
ellipse(310,201,4.7,4.7);
fill(#000000);
ellipse(285,190,7,7);

fill(46,134,126);
rect(200,325,100,5.3);
