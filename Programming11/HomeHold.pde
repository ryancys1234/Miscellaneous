int horizontal = 1;
color darkbluepink,lightbluepink;

void setup() {
  size(1280,720);
  
  darkbluepink = color(#040043);
  lightbluepink = color(#8067E8);
  noLoop();
}

void draw() {
  gradient(0,0,1280,720,darkbluepink,lightbluepink,horizontal);
  
  fill(#1A105F);
  noStroke();
  triangle(0,600,200,450,340,580);
  
  fill(#FFFFFF);
  ellipse(492,182,3,3);
  ellipse(257,293,3,3);
}

void gradient(int x, int y, float w, float h, color darkbluepink, color lightbluepink, int axis) {
  noFill();
  
  if (axis == horizontal) {
    for (int i = y;i <= y+h;i++) {
      float inter = map(i,y,y+h,0,1);
      color c = lerpColor(darkbluepink,lightbluepink,inter);
      stroke(c);
      line(x,i,x+w,i);
    }
  }
}

