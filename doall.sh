## ROSEN
BASE=rosen8-dog
echo -e "4\n1\n1" | ./rosen.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp delta.png ../RESOURCES/${BASE}-d-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n1" | ./rosen.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp delta.png ../RESOURCES/${BASE}-d-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n1" | ./rosen.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp delta.png ../RESOURCES/${BASE}-d-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

BASE=rosen8-line
echo -e "4\n1\n2" | ./rosen.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp backtracks.png ../RESOURCES/${BASE}-b-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n2" | ./rosen.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp backtracks.png ../RESOURCES/${BASE}-b-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n2" | ./rosen.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp backtracks.png ../RESOURCES/${BASE}-b-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

echo -e "4\n4\n" | ./rosen.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nohess.png
cp backtracks.png ../RESOURCES/${BASE}-b-nohess.png
cp answer.out ../RESOURCES/${BASE}-nohess.out

echo -e "4\n5\n" | ./rosen.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nograd.png
cp backtracks.png ../RESOURCES/${BASE}-b-nograd.png
cp answer.out ../RESOURCES/${BASE}-nograd.out

## POWELL

BASE=powell8-dog
echo -e "4\n1\n1" | ./powell.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp delta.png ../RESOURCES/${BASE}-d-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n1" | ./powell.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp delta.png ../RESOURCES/${BASE}-d-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n1" | ./powell.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp delta.png ../RESOURCES/${BASE}-d-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

BASE=powell8-line
echo -e "4\n1\n2" | ./powell.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp backtracks.png ../RESOURCES/${BASE}-b-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n2" | ./powell.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp backtracks.png ../RESOURCES/${BASE}-b-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n2" | ./powell.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp backtracks.png ../RESOURCES/${BASE}-b-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

echo -e "4\n4\n" | ./powell.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nohess.png
cp backtracks.png ../RESOURCES/${BASE}-b-nohess.png
cp answer.out ../RESOURCES/${BASE}-nohess.out

echo -e "4\n5\n" | ./powell.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nograd.png
cp backtracks.png ../RESOURCES/${BASE}-b-nograd.png
cp answer.out ../RESOURCES/${BASE}-nograd.out

## WOOD

BASE=wood8-dog
echo -e "4\n1\n1" | ./wood.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp delta.png ../RESOURCES/${BASE}-d-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n1" | ./wood.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp delta.png ../RESOURCES/${BASE}-d-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n1" | ./wood.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp delta.png ../RESOURCES/${BASE}-d-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

BASE=wood8-line
echo -e "4\n1\n2" | ./wood.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp backtracks.png ../RESOURCES/${BASE}-b-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n2" | ./wood.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp backtracks.png ../RESOURCES/${BASE}-b-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n2" | ./wood.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp backtracks.png ../RESOURCES/${BASE}-b-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

echo -e "4\n4\n" | ./wood.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nohess.png
cp backtracks.png ../RESOURCES/${BASE}-b-nohess.png
cp answer.out ../RESOURCES/${BASE}-nohess.out

echo -e "4\n5\n" | ./wood.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nograd.png
cp backtracks.png ../RESOURCES/${BASE}-b-nograd.png
cp answer.out ../RESOURCES/${BASE}-nograd.out

### quad!

## ROSEN
BASE=rosen16-dog
echo -e "4\n1\n1" | ./rosen.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp delta.png ../RESOURCES/${BASE}-d-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n1" | ./rosen.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp delta.png ../RESOURCES/${BASE}-d-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n1" | ./rosen.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp delta.png ../RESOURCES/${BASE}-d-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

BASE=rosen16-line
echo -e "4\n1\n2" | ./rosen.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp backtracks.png ../RESOURCES/${BASE}-b-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n2" | ./rosen.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp backtracks.png ../RESOURCES/${BASE}-b-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n2" | ./rosen.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp backtracks.png ../RESOURCES/${BASE}-b-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

echo -e "4\n4\n" | ./rosen.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nohess.png
cp backtracks.png ../RESOURCES/${BASE}-b-nohess.png
cp answer.out ../RESOURCES/${BASE}-nohess.out

echo -e "4\n5\n" | ./rosen.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nograd.png
cp backtracks.png ../RESOURCES/${BASE}-b-nograd.png
cp answer.out ../RESOURCES/${BASE}-nograd.out

## POWELL

BASE=powell16-dog
echo -e "4\n1\n1" | ./powell.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp delta.png ../RESOURCES/${BASE}-d-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n1" | ./powell.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp delta.png ../RESOURCES/${BASE}-d-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n1" | ./powell.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp delta.png ../RESOURCES/${BASE}-d-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

BASE=powell16-line
echo -e "4\n1\n2" | ./powell.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp backtracks.png ../RESOURCES/${BASE}-b-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n2" | ./powell.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp backtracks.png ../RESOURCES/${BASE}-b-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n2" | ./powell.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp backtracks.png ../RESOURCES/${BASE}-b-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

echo -e "4\n4\n" | ./powell.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nohess.png
cp backtracks.png ../RESOURCES/${BASE}-b-nohess.png
cp answer.out ../RESOURCES/${BASE}-nohess.out

echo -e "4\n5\n" | ./powell.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nograd.png
cp backtracks.png ../RESOURCES/${BASE}-b-nograd.png
cp answer.out ../RESOURCES/${BASE}-nograd.out

## WOOD

BASE=wood16-dog
echo -e "4\n1\n1" | ./wood.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp delta.png ../RESOURCES/${BASE}-d-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n1" | ./wood.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp delta.png ../RESOURCES/${BASE}-d-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n1" | ./wood.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp delta.png ../RESOURCES/${BASE}-d-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

BASE=wood16-line
echo -e "4\n1\n2" | ./wood.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-ana.png
cp backtracks.png ../RESOURCES/${BASE}-b-ana.png
cp answer.out ../RESOURCES/${BASE}-ana.out
echo -e "4\n2\n2" | ./wood.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-grad.png
cp backtracks.png ../RESOURCES/${BASE}-b-grad.png
cp answer.out ../RESOURCES/${BASE}-grad.out
echo -e "4\n3\n2" | ./wood.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-fake.png
cp backtracks.png ../RESOURCES/${BASE}-b-fake.png
cp answer.out ../RESOURCES/${BASE}-fake.out

echo -e "4\n4\n" | ./wood.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nohess.png
cp backtracks.png ../RESOURCES/${BASE}-b-nohess.png
cp answer.out ../RESOURCES/${BASE}-nohess.out

echo -e "4\n5\n" | ./wood.16.exe  	
gnuplot makeplots.plot
cp function.png ../RESOURCES/${BASE}-f-nograd.png
cp backtracks.png ../RESOURCES/${BASE}-b-nograd.png
cp answer.out ../RESOURCES/${BASE}-nograd.out
