using System;


class WolvesGWO{
	private static Random random = new Random();

	private const int AGENT = 100;
	private const int DIMENSION = 2;
	private const int MAX_IT = 50000;

	private static int alpha;
	private static int beta;
	private static int delta;
	private static double[,]wolves = new double[AGENT,DIMENSION];
	private static double[]alphaPos = new double[DIMENSION];
	private static double[]betaPos = new double[DIMENSION];
	private static double[]deltaPos = new double[DIMENSION];

	private static double getRandom(double min, double max){
		return random.NextDouble() * (max - min) + min;
	}

	public static double[] optimize(){
		initialize();
		for(int i=0;i<WolvesGWO.MAX_IT;i++){
			update_wolves();
			update(i);
		}
		Console.WriteLine("f(x,y) = " + function(WolvesGWO.alphaPos) + " || x = " + WolvesGWO.alphaPos[0] + " || y = " + WolvesGWO.alphaPos[1]);
		return WolvesGWO.alphaPos;
	}

	private static void initialize(){
		for(int i=0;i<AGENT;i++){
			WolvesGWO.wolves[i,0] = getRandom(-10,10);
			WolvesGWO.wolves[i,1] = getRandom(-10,10);
		}

		WolvesGWO.alpha = elitism(Double.MaxValue);
		WolvesGWO.beta = elitism(WolvesGWO.alpha);
		WolvesGWO.delta = elitism(WolvesGWO.beta);

	}

	private static int elitism(double found_maximum){
		double[] WolvesGWOVector = new double[2];

		double maximum = fitness(WolvesGWOVector);
		int index_maximum = 0;
			for(int i=1;i<AGENT;i++){
				  Buffer.BlockCopy(WolvesGWO.wolves,i*8*WolvesGWO.DIMENSION,WolvesGWOVector,0,WolvesGWO.DIMENSION*8);
					if(fitness(WolvesGWOVector) > maximum && fitness(WolvesGWOVector) < found_maximum){
							index_maximum = i;
							maximum = fitness(WolvesGWOVector);
						}
			 }
	 	return index_maximum;
	}

	private static void update_wolves(){
			double result;
			int aux_index_alpha;
			int aux_index_beta;
			int aux_index_delta;
			double[] WolvesGWOVector = new double[2];

			for(int i=0;i<WolvesGWO.AGENT;i++){
					Buffer.BlockCopy(WolvesGWO.wolves,i*8*WolvesGWO.DIMENSION,WolvesGWOVector,0,WolvesGWO.DIMENSION*8);
					result = fitness(WolvesGWOVector);
					if(result > fitness(WolvesGWO.alphaPos)){
								aux_index_alpha = WolvesGWO.alpha;
								aux_index_beta = WolvesGWO.beta;

								Buffer.BlockCopy(WolvesGWO.wolves,i*8*WolvesGWO.DIMENSION,WolvesGWO.alphaPos,0,WolvesGWO.DIMENSION*8);
								Buffer.BlockCopy(WolvesGWO.wolves,WolvesGWO.alpha*8*WolvesGWO.DIMENSION,WolvesGWO.betaPos	,0,WolvesGWO.DIMENSION*8);
								Buffer.BlockCopy(WolvesGWO.wolves,WolvesGWO.beta*8*WolvesGWO.DIMENSION,WolvesGWO.deltaPos,0,WolvesGWO.DIMENSION*8);

								WolvesGWO.alpha = i;
                                                                WolvesGWO.beta  = aux_index_alpha;
                                                                WolvesGWO.delta = aux_index_beta;

						

					}
					else if(result > fitness(WolvesGWO.betaPos) && result < fitness(WolvesGWO.alphaPos)){
								aux_index_beta = WolvesGWO.beta;
								aux_index_delta = WolvesGWO.delta;

								Buffer.BlockCopy(WolvesGWO.wolves,i*8*WolvesGWO.DIMENSION,WolvesGWO.betaPos,0,WolvesGWO.DIMENSION*8);
								Buffer.BlockCopy(WolvesGWO.wolves,WolvesGWO.beta*8*WolvesGWO.DIMENSION,WolvesGWO.deltaPos,0,WolvesGWO.DIMENSION*8);

								WolvesGWO.beta = i;
                                                                WolvesGWO.delta = aux_index_beta;
					}
					else if (result > fitness(WolvesGWO.deltaPos) && result < fitness(WolvesGWO.betaPos)){
								WolvesGWO.delta = i;
								Buffer.BlockCopy(WolvesGWO.wolves,i*8*WolvesGWO.DIMENSION,WolvesGWO.deltaPos,0,WolvesGWO.DIMENSION*8);

					}
			}

	}

	private static int update_a(int iteration){
			return 2 - (Convert.ToDouble(iteration)*(2/( Convert.ToDouble(WolvesGWO.MAX_IT)-1)));
	}

	private static void update(int iteration){
			double a = update_a(iteration);
			for(int i=0;i<WolvesGWO.AGENT;i++){
						for(int j=0;j<DIMENSION;j++){

							double r1 = random.NextDouble();
							double r2 = random.NextDouble();
							double a1 = 2*a*r1 - a;
							double c1 = 2*r2;
							double dAlpha = (c1*WolvesGWO.alphaPos[j] - WolvesGWO.wolves[i,j]);
							double x1 = (WolvesGWO.alphaPos[j] - (a1*dAlpha));

							r1 = random.NextDouble();
							r2 = random.NextDouble();
							double a2 = 2*a*r1 - a;
							double c2 = 2*r2;
							double dBeta = (c2*WolvesGWO.betaPos[j] - WolvesGWO.wolves[i,j]);
							double x2 = (WolvesGWO.betaPos[j] - (a2*dBeta));

							r1 = random.NextDouble();
							r2 = random.NextDouble();
							double a3 = 2*a*r1 - a;
							double c3 = 2*r2;
							double dDelta = (c3*WolvesGWO.deltaPos[j] - WolvesGWO.wolves[i,j]);
							double x3 = (WolvesGWO.deltaPos[j] - (a3*dDelta));

							WolvesGWO.wolves[i,j] = (x1+x2+x3)/3;

					}

			}


	}

	private static double fitness(double[] x){
		if(function(x) >= 0)
				return 1/(1+function(x));
		else
				return 1 + (-1*(function(x)));
	}

	private static double function(double[] x){
		//return 0.26*(Math.Pow(x[0],2) + Math.Pow(x[1],2)) - 0.48*x[0]*x[1];  // Matyas Function
		//return Math.Pow(x[0] + 10,2) + Math.Pow(x[1] + 10,2) + Math.Exp(-Math.Pow(x[0],2)-Math.Pow(x[1],2)); // Brent  Function
		//return 100*(Math.Sqrt(Math.Abs(x[1] - 0.01*Math.Pow(x[0],2)))) + 0.01*Math.Abs(x[0] + 10); // Bukin N.6 Function
		return Math.Pow(Math.Pow(x[0],2) + x[1]- 11,2) + Math.Pow(x[0] + Math.Pow(x[1],2) - 7,2);  // Himmelblau Function
	}




}
