int n = 6;
float scale = 0.9;
ArrayList<Robot> robots;
boolean convexHull = false;
boolean random = false;
boolean pointing = false;
float dia;

void setup()
{
  size(800,800);
  frameRate(30);
  dia = width * scale;
  robots = new ArrayList<Robot>();
  
  if (random) randomDistribution();
  else evenDistribution();
  
 
 
}

void draw()
{
  background(255);
  fill(255);
  circle(height/2,width/2,dia);
  circle(mouseX,mouseY,5);
  for (int i=0; i<n; i++){
    fill(0);
    circle(robots.get(i).pos.x,robots.get(i).pos.y, 5);
    //robots[i].update();
    robots.get(i).drawLine();
  
  }
  if (convexHull) drawConvexHull();
  if (pointing) drawPointing();
  //println(n);
  //save("test2.jpg");
}


class Robot {
  PVector pos = new PVector(0,0);
  float conc = 0.0;
  PVector grad = new PVector(0,0);
  
  Robot(float x, float y) {
    pos = new PVector(x,y);
  }
  
  void update(){
    
  }
  
  void drawLine(){
    //line(pos.x,pos.y, mouseX, mouseY);
    //float slope = (mouseY - pos.y) / (mouseX - pos.x);
    //float slope = (pos.y-mouseY) / (pos.x - mouseX);
    //slope = -1 / slope;
    //float b = pos.y/(slope * pos.x);
    //float left_y = (slope * 0) + b;
    //float right_y = (slope * width) + b;
    //line(pos.x,pos.y, 0, left_y);
    //line(pos.x,pos.y, width, right_y);
    PVector mouse = new PVector(mouseX,mouseY);
    PVector vec = PVector.sub(mouse, pos);
    vec.setMag(width);
    PVector right = vec.copy();
    PVector left = vec.copy();
    right.rotate(radians(90));
    left.rotate(radians(-90));
    right.add(pos);
    left.add(pos);
    line(pos.x,pos.y,right.x,right.y);
    line(pos.x,pos.y,left.x,left.y);
    
    
  }
  
 
  
}

void keyPressed() {
   if (key=='c') convexHull = !convexHull;
   if (key=='=') {
     n++;
     setup();
   }
   if (key=='-'){
     n--;
     if (n<1) n = 1;
     setup();
   }
   if (key=='r'){
     random = true;
     setup();
   }
   if (key=='e'){
     random = false;
     setup();
   }
   if (key=='l'){
     pointing = !pointing;
     setup();
   }
}

void drawConvexHull(){
  for (int i=0; i<n; i++){
    int other;
      if (i==0) {  
        other = n-1;
      }
      else{
        other = i-1; 
      }
      PVector otherPos = robots.get(other).pos;
      line(robots.get(i).pos.x,robots.get(i).pos.y,otherPos.x,otherPos.y);
  }
}

void drawPointing(){
  for (int i=0; i<n; i++){
   line(robots.get(i).pos.x,robots.get(i).pos.y, mouseX, mouseY);
  }
}

void evenDistribution(){
  float angle_delta = radians(360.0/n);
  float r = ((width * scale) / 2) / 2;
 
  if (n > 1){
    for (int i=0; i<n; i++){
      float x = r * cos(i*angle_delta) + width/2;
      float y = r * sin(i*angle_delta) + width/2;
      robots.add(new Robot(x,y));
    }
  }
  else robots.add(new Robot(width/2, height/2));

  
}

void randomDistribution(){
  robots.add(new Robot(width/2,height/2));
  for (int i=0; i<n-1; i++){
    boolean done = false;
    float x = 0;
    float y = 0;
    float r = dia / 2;
    while(!done){
      x = random(0, width);
      y = random(0, height);
      if (dist(x,y,width/2,height/2) <= r) done = true;
    }
    robots.add(new Robot(x,y));
  }
}
