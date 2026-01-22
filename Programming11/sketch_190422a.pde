int red=255;
int green=255;
int blue=255;

void setup() {
  size(800,600);
}

void draw() {
  background(red,green,blue);
 
  noStroke();
  if (dist(400,150,mouseX,mouseY) <= 50) {
    fill(#2725CE);
  } else {
    fill(#5D5BF5);
  }
  
  ellipse(400,150,100,100);
 
  noStroke();
  if (dist(400,300,mouseX,mouseY) <= 50) {
    fill(255,0,0);
  } else {
    fill(#E55858);
  }
  
  ellipse(400,300,100,100);
  
  noStroke();
  if (dist(400,450,mouseX,mouseY) <= 50) {
    fill(#0FD33B);
  } else {
    fill(#48D18F);
  }
  
  ellipse(400,450,100,100);
}

void mouseReleased() {
  if (dist(400,150,mouseX,mouseY) <=50) {
    red=0;
    green=0;
    blue=255;
  }
}
