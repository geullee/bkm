package kr.geul.bkm;

import kr.geul.options.exception.DuplicateOptionsException;
import kr.geul.options.exception.InconsistentOptionException;
import kr.geul.options.option.Option;
import kr.geul.options.structure.OptionCurve;

public class BKMEstimator {

	double underlyingPrice, timeToMaturity, riskFreeRate, dividendRate,
	v, w, x, mu, vol, skew, kurt,
	extrapolationLeftEnd, extrapolationRightEnd, strikePriceGap;
	OptionCurve callCurve, putCurve;

	public void setOptions(OptionCurve curve) throws DuplicateOptionsException, 
	InconsistentOptionException {

		double[] variables = curve.getVariableArray();
		underlyingPrice = variables[0];
		timeToMaturity = variables[1];
		riskFreeRate = variables[2];
		dividendRate = variables[3];
		callCurve = curve.getCallCurve();
		putCurve = curve.getPutCurve();
		
		System.out
				.println("OptionCurve passed to BKMEstimator: "
						+ putCurve.size() + " puts, " + callCurve.size()
						+ " calls"); 
		
		System.out
				.println("Kmin_put: "
						+ putCurve.getStrikePrices()[0]
						+ ", Kmax_put: "
						+ putCurve.getStrikePrices()[putCurve.getStrikePrices().length - 1]
						+ ", Kmin_call: "
						+ callCurve.getStrikePrices()[0]
						+ ", Kmax_call: "
						+ callCurve.getStrikePrices()[callCurve.getStrikePrices().length - 1]); 
		
	}
	
	public double[] getEstimates() {

		getVWX();
		getMoments();

		double[] estimates = {vol, skew, kurt, v, w, x, mu};
		return estimates;

	}

	private void getMoments() {

		double ert = Math.exp(riskFreeRate * timeToMaturity);
		mu = ert - 1 - (ert * v / 2) - (ert * w / 6) - (ert * x / 24);	
		vol = Math.sqrt(((ert * v) - Math.pow(mu, 2)) / timeToMaturity);
		skew = ((ert * w) - (3 * mu * ert * v) + (2 * Math.pow(mu, 3))) / 
				Math.pow((ert * v) - Math.pow(mu, 2), 1.5); 
		kurt = ((ert * x) - (4 * mu * ert * w) + (6 * ert * Math.pow(mu, 2) * v) - 
				(3 * Math.pow(mu, 4))) / Math.pow((ert * v) - Math.pow(mu, 2), 2); 

	}

	private double getTrapezoidalIntegrationValue(double[] xValues, double[] yValues) {

		double totalArea = 0.0;

		for (int i = 1; i < xValues.length; i++) {	

			double area = (xValues[i] - xValues[i - 1]) * 
					((yValues[i] + yValues[i - 1]) * 0.5); 
			totalArea += area;

		}

		return totalArea;

	}

	private void getVWX() {
		
		double[] vWeightedPrices_put = new double[putCurve.size()],
				wWeightedPrices_put = new double[putCurve.size()],
				xWeightedPrices_put = new double[putCurve.size()],
				vWeightedPrices_call = new double[callCurve.size()],
				wWeightedPrices_call = new double[callCurve.size()],
				xWeightedPrices_call = new double[callCurve.size()];

		for (int i = 0; i < putCurve.size(); i++) {

			Option option = putCurve.get(i);

			double optionPrice = option.getOptionPrice(),
					strikePrice = option.getStrikePrice(),
					vWeightedPrice = (2.0 * (1.0 - Math.log(strikePrice / (underlyingPrice * 
							Math.exp(-dividendRate * timeToMaturity))))) / 
					(Math.pow(strikePrice, 2.0)) * optionPrice,
					wWeightedPrice = ((6.0 * Math.log(strikePrice / (underlyingPrice * 
							Math.exp(-dividendRate * timeToMaturity)))) - 
							(3.0 * Math.pow(Math.log(strikePrice / (underlyingPrice * 
									Math.exp(-dividendRate * timeToMaturity))), 2.0))) / 
							(Math.pow(strikePrice, 2.0)) * optionPrice,
					xWeightedPrice = ((12.0 * Math.pow(Math.log(strikePrice / (underlyingPrice * 
							Math.exp(-dividendRate * timeToMaturity))), 2.0)) - 
							(4.0 * Math.pow(Math.log(strikePrice / (underlyingPrice * 
									Math.exp(-dividendRate * timeToMaturity))), 3.0))) 
							/ (Math.pow(strikePrice, 2.0)) * optionPrice;
			
			vWeightedPrices_put[i] = vWeightedPrice;
			wWeightedPrices_put[i] = wWeightedPrice;
			xWeightedPrices_put[i] = xWeightedPrice;

		}
		
		for (int i = 0; i < callCurve.size(); i++) {

			Option option = callCurve.get(i);
			
			double optionPrice = option.getOptionPrice(),
					strikePrice = option.getStrikePrice(),
					vWeightedPrice = (2.0 * (1.0 - Math.log(strikePrice / (underlyingPrice * 
							Math.exp(-dividendRate * timeToMaturity))))) / 
					(Math.pow(strikePrice, 2.0)) * optionPrice,
					wWeightedPrice = ((6.0 * Math.log(strikePrice / (underlyingPrice * 
							Math.exp(-dividendRate * timeToMaturity)))) - 
							(3.0 * Math.pow(Math.log(strikePrice / (underlyingPrice * 
									Math.exp(-dividendRate * timeToMaturity))), 2.0))) / 
							(Math.pow(strikePrice, 2.0)) * optionPrice,
							xWeightedPrice = ((12.0 * Math.pow(Math.log(strikePrice / (underlyingPrice * 
									Math.exp(-dividendRate * timeToMaturity))), 2.0)) - 
									(4.0 * Math.pow(Math.log(strikePrice / (underlyingPrice * 
											Math.exp(-dividendRate * timeToMaturity))), 3.0))) 
									/ (Math.pow(strikePrice, 2.0)) * optionPrice;

			vWeightedPrices_call[i] = vWeightedPrice;
			wWeightedPrices_call[i] = wWeightedPrice;
			xWeightedPrices_call[i] = xWeightedPrice;

		}

		v = getTrapezoidalIntegrationValue(putCurve.getStrikePrices(), vWeightedPrices_put)
				+ getTrapezoidalIntegrationValue(callCurve.getStrikePrices(), vWeightedPrices_call);
		w = getTrapezoidalIntegrationValue(putCurve.getStrikePrices(), wWeightedPrices_put)
				+ getTrapezoidalIntegrationValue(callCurve.getStrikePrices(), wWeightedPrices_call);
		x = getTrapezoidalIntegrationValue(putCurve.getStrikePrices(), xWeightedPrices_put)
				+ getTrapezoidalIntegrationValue(callCurve.getStrikePrices(), xWeightedPrices_call);

	}

}
