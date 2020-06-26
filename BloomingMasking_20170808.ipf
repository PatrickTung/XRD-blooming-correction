#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Peak AutoFind>

//This was my attempt at writing a procedure to find peaks along a 1D wave
//Didn't seem to work because sometimes the PeakWidth would be NaN
//This affected how I recalculated rangeStart through the loop, and this would result in a infinite loop
function FindPeaks(input, minLevel, boxSize)
	wave input	//input image
	variable minLevel	//the minimum intensity to even be considered a peak
	variable boxSize	//sets box size for sliding average. Search FindPeak command for more info. Test how this affects peak finding with test data. 30 works.
	
	make/O/N=(dimsize(input,0)) peakWaveLine	//make a wave to hold a row of data from the 2D image
	duplicate/O input, peakWave	//make a wave that shows where the peaks occur
	peakWave = 0	//may be faster if multithreaded
	
	variable i,j,rangeStart
	for(i=0;i<dimsize(input,1);i+=1)	//loop through the columns of the 2D image
		peakWaveLine = input[p][i]	//extract a row of data from 2D image
		rangeStart = 0	//reset the value of rangeStart for the next line
		do
			FindPeak /Q /M=(minLevel) /B=(boxSize) /R=[rangeStart,dimsize(peakWaveLine,0)] peakWaveLine	//find the peak within the range specified
			if(numtype(V_PeakLoc)==0)
				peakWave[round(V_PeakLoc)][i] = V_PeakVal		//set the location of the peak as the intensity of the peak
				rangeStart = V_PeakLoc + V_PeakWidth 	//change the value of the start of the peak search
				if(numtype(rangeStart) ==2)
					print i
				endif
			endif
		while(V_flag == 0)
		
	endfor
		
	return 0	
end

//Function will find peaks along a 1D wave with the use of in-built functions, and return a 2D wave with peaks found in one direction with a value of 1
//To include in-built functions, type #include <Peak AutoFind> at the top of the procedure file
Function FindPeaksPro(inWave, minLevel, maxPeaks,minPeakPercent)
	wave inWave		//input wave for 2D image
	Variable minLevel		//minimum intensity accepted for a peak
	Variable maxPeaks	//maximum number of peaks expected within the image
	Variable minPeakPercent	//the minimum intensity accepted as a % of max intensity
	
	make/O/N=(dimsize(inWave,0)) peakWaveLine	//make a wave to hold a row of data from the 2D image
	duplicate/O inWave, peakWave	//make a wave that shows where the peaks occur
	peakWave = 0	//may be faster if multithreaded
	
	//Set up all the variables required for using the built-in peak finding function
	Variable pBegin=0, pEnd= dimsize(inWave,0)-1
	Variable/C estimates
	Variable noiselevel, smoothingFactor
	Wave WA_PeakCentersX, WA_PeakCentersY		//reference for wave that holds peak center coordinates
	make/O/N=(dimsize(inwave,0)) wx		//make a wave for the x-coordinates of the wave
	wx = p	//just make the values of the wave the same as the row number
	
	variable i,j
	for(i=0;i<dimsize(inWave,1);i+=1)	//loop through the columns of the 2D image
		peakWaveLine = inWave[p][i]
		wavestats/Q peakWaveLine
		if(V_max > minLevel)
			estimates= EstPeakNoiseAndSmfact(peakWaveLine,pBegin, pEnd)
			noiselevel=real(estimates)
			smoothingFactor=imag(estimates)
			AutoFindPeaksWorker(peakWaveLine, wx, pBegin, pEnd, maxPeaks, minPeakPercent, noiseLevel, smoothingFactor)
//			Doupdate		//uncomment this to see diagnostics in real-time
			for(j=0;j<dimsize(WA_PeakCentersX,0);j+=1)		//loop through the rows of the resulting peak-finding waves
				variable peakCent = round(WA_PeakCentersX[j])
//				peakWave[peakCent][i] = WA_PeakCentersY[j]		//set the location of the peak as the intensity of the peak
				peakWave[peakCent][i] += 1	//set the location of the peak as 1
			endfor
		endif
		if(mod(i,250)==0)
			print i
		endif
	endfor
	
	////////////////////////////////////////////////////////////////Vertical direction///////////////////////////////////////////////////////////////////
	
	redimension/N=(dimsize(inWave,1)) peakWaveLine
	
	for(i=0;i<dimsize(inWave,0);i+=1)	//loop through the columns of the 2D image
		peakWaveLine = inWave[i][p]
		wavestats/Q peakWaveLine
		if(V_max > minLevel)
			estimates= EstPeakNoiseAndSmfact(peakWaveLine,pBegin, pEnd)
			noiselevel=real(estimates)
			smoothingFactor=imag(estimates)
			AutoFindPeaksWorker(peakWaveLine, wx, pBegin, pEnd, maxPeaks, minPeakPercent, noiseLevel, smoothingFactor)
//			Doupdate		//uncomment this to see diagnostics in real-time
			for(j=0;j<dimsize(WA_PeakCentersX,0);j+=1)		//loop through the rows of the resulting peak-finding waves
				peakCent = round(WA_PeakCentersX[j])
//				peakWave[i][peakCent] += WA_PeakCentersY[j]		//set the location of the peak as the intensity of the peak
				peakWave[i][peakCent] += 1	//set the location of the peak as 1
			endfor
		endif
		if(mod(i,250)==0)
			print i
		endif
	endfor
	
	killwaves wx, WA_PeakCentersX, WA_PeakCentersY
	
	return 0
	
end

function StreakDetectionSobel(inWave, kernel, maxThresh, outWave)
	wave inWave
	wave kernel
	variable maxThresh	// The maximum intensity that will be used. e.g. 1250 for ESRF 2013 data
	string outWave
	
	duplicate/O inWave, $outWave, maskWave
	wave OW = $outWave
	OW = inWave[p][q] > maxThresh ? maxThresh : OW[p][q]		//Make a histogram of the image and with the Histogram function and then decide the threshold value	
	maskWave = inWave > maxThresh ? 1 : 0
	duplicate/FREE inWave, OWx, OWy		//waves to hold the values from horizontal and vertical sobel operations
	
	//Magic happens here: image is greyscaled, blur applied, then kernel convolution
	ImageTransform/O convert2gray OW	//grayscale the image and overwrite
	MatrixFilter gauss OW					//apply a gaussian blur
//	MatrixFilter avg OW					//apply a mean blur
	MatrixConvolve kernel, OW				//apply matrix convolution with the specified kernel

//	Sobelwave = atan(SobelWavex/sobelWavey)*180/pi

	redimension/S OW					//redimension to single point floating number, so that NaN's can be used
//	OW = inWave > maxThresh ? Nan : OW
end


function StreakDetectionVert(OGWave,inWave,maxThresh, minThresh)
	wave OGWave		//original detector image
	wave inWave			//input wave
	variable maxThresh	// The maximum intensity that will be used. e.g. 1250 for ESRF 2013 data
	variable minThresh	//the minimum amount of consecutive pixels to be considered a streak of abnormality
	
	duplicate/O inWave, mask		//make a mask wave that will be the mask for the 2D image
	multithread mask = 1
	
	////////////////////Looping through columns  of rows///////////////////////
	variable i,j,k, counter
	make/O/N=(dimsize(mask,1)) maskLine, maskTracker = 0	//make a wave to hold the data that will be analysed

	variable maskTrack
	for(i=0;i<dimsize(mask,0);i+=1)		//loop through the rows
		maskLine = inWave[i][p]	//extract a column of data
	
		for(j=1;j<dimsize(mask,1);j+=1)	//compare each point to its neigbouring points from bottom to top
			if(OGWave[i][j] < maxThresh && OGWave[i][j-1] >= maxThresh)		//if the pixel of interest is connected to a high intensity
				do		//search the rest of the line to see if there are consecutive pixels of the same value
					if(maskLine[j] == maskLine[j-1] && OGWave[i][j] < maxThresh)	//if the previous pixel is the same AND the pixel from the original wave is not a high intensity
						maskTrack += 1	//increment the tracker for the pixel along the row by 1
						if(maskTrack > minThresh && maskLine[j] == 0)	//if the amount of consecutive pixels is over the minimum threshold
//							mask[i][(j-minThresh), j] = 0		//make that section equal to 0
							mask[i][] = 0	//mask the whole line
						elseif(maskTrack > minThresh && maskLine[j] == 255)
							mask[i][] = 255	//mask the whole line
						endif
						j += 1
					else
						maskTrack = 0		//reset the counter
					endif
				while(maskTrack > 0 && j <dimsize(mask,1))
			endif
		endfor
		
		maskTrack = 0		//reset the counter for the next line
	endfor

	for(i=0;i<dimsize(mask,0);i+=1)		//loop through the rows
		maskLine = inWave[i][p]
		for(j=dimsize(mask,1)-2; j>=0; j-=1)	//compare each point to its neigbouring points from top to bottom
			if(OGWave[i][j] < maxThresh && OGWave[i][j+1] >= maxThresh)	
				do		//search the rest of the line to see if there are consecutive pixels of the same value
					if(maskLine[j] == maskLine[j+1] && OGWave[i][j] < maxThresh)	//if the previous pixel is the same AND the pixel from the original wave is not a high intensity
						maskTrack += 1	//increment the tracker for the pixel along the row by 1
						if(maskTrack > minThresh && maskLine[j] == 0)	//if the amount of consecutive pixels is over the minimum threshold
//							mask[i][(j-minThresh), j] = 0		//make that section equal to 0
							mask[i][] = 0	//mask the whole line
						elseif(maskTrack > minThresh && maskLine[j] == 255)
							mask[i][] = 255	//mask the whole line
						endif
						j -= 1
					else
						maskTrack = 0		//reset the counter
					endif
				while(maskTrack > 0 && j >= 0)
			endif
		endfor

		maskTrack = 0		//reset the counter for the next line
	endfor
	
	//fill in the empty space between the lines to make a solid line
	for(i=0;i<dimsize(mask,0);i+=1)
		if(mask[i][0] == 255)
			mask[i-1][] = 0
			mask[i-2][] = 0
		endif
	endfor
	
	mask = mask == 255 ? 0 : mask
	//killwaves maskTracker, maskLine
	
	return 0
end


function StreakDetectionHori(inWave,minThresh)
	wave inWave			//input wave
	variable minThresh	//the minimum amount of consecutive pixels to be considered a streak of abnormality
	
	duplicate/O inWave, mask		//make a mask wave that will be the mask for the 2D image
	multithread mask = 1
	
	////////////////////Looping through rows of columns///////////////////////
	variable i,j,k
	make/O/N=(dimsize(mask,0)) maskLine, maskTracker = 0	//make a wave to hold the data that will be analysed

	variable maskTrack
	for(i=0;i<dimsize(mask,1);i+=1)		//loop through the rows of columns
		maskLine = inWave[p][i]	//extract a row of data
		for(j=1;j<dimsize(mask,0);j+=1)	//compare each point to its neigbouring points. ignore the end points
			if(maskLine[j] == maskLine[j-1])	
				maskTrack += 1	//increment the tracker for the pixel along the row by 1
				if(maskTrack > minThresh)	//if the amount of consecutive pixels is over the minimum threshold
					for(k=j;k>j-minThresh;k-=1)	//loop backwards from the current pixel we're looking at, up until the amount where the consecutive counts started
						mask[k][i] = 0
					endfor
				endif
			else
				maskTrack = 0		//reset the counter
			endif
		endfor
		maskTrack = 0
	endfor

	//killwaves maskTracker, maskLine
	
	return 0
end

//Makes a mask that can be used to remove streaks that run along columns of data
function StreakMask(inWave, maxThresh, minStreakLength, peakMask, maskName)
	wave inWave				//2D detector image
	variable maxThresh		//maximum intensity that you are willing to accept as background diffuse scattering e.g. 1250 counts, anything over will be masked
	variable minStreakLength		//the minimum amount of consecutive pixels to be considered a streak of abnormality
	wave peakMask
	string maskName
	
	//make the sobel kernel that finds edges along the columns of the detector image
	make/FREE/N=25 sobelKernelVert5x5
	sobelKernelVert5x5 = {-2,-1,0,1,2, -3,-2,0,2,3, -4,-3,0,3,4, -3,-2,0,2,3, -2,-1,0,1,2}	//put all in one column
	redimension/n=(5,5) sobelKernelVert5x5		//then redimension to 5 columns to make a 5x5 matrix
	
	//Magic happens here: image is greyscaled, blur applied, then kernel convolution
	duplicate/FREE inWave, sobelWave
	multithread sobelWave = inWave > maxThresh ? maxThresh : sobelWave		//if pixel are greater than the max threshold, then make those pixels in sobelWave equal to max threshold
	ImageTransform/O convert2gray sobelWave	//grayscale the image and overwrite
	MatrixFilter gauss sobelWave					//apply a gaussian blur
//	MatrixFilter avg OW							//apply a mean blur
	MatrixConvolve sobelKernelVert5x5, sobelWave	//apply matrix convolution with the specified kernel
//	Sobelwave = atan(SobelWavex/sobelWavey)*180/pi	//to find the orientation of the edge
	
	duplicate/O inWave, $maskName		//make a mask wave that will be the mask for the 2D image
	wave mask = $maskName
	redimension/S mask	//redimension to single floating point so that it can have NaN's
	multithread mask = 1
	
	////////////////////Looping through columns  of rows///////////////////////
	variable i,j,k, counter
	make/FREE/N=(dimsize(mask,1)) maskLine	//make a wave to hold the data that will be analysed

	variable maskTrack
	for(i=0;i<dimsize(mask,0);i+=1)		//loop through the rows
		maskLine = sobelWave[i][p]	//extract a column of data
	
		for(j=1;j<dimsize(mask,1);j+=1)	//compare each point to its neigbouring points from bottom to top
//			if(inWave[i][j] < maxThresh && inWave[i][j-1] >= maxThresh)		//if the pixel of interest is connected to a high intensity
			if(inWave[i][j] < maxThresh && peakMask[i][j-1] == 0)		//if the pixel of interest is connected to a high intensity
				do		//search the rest of the line to see if there are consecutive pixels of the same value
					if(maskLine[j] == maskLine[j-1] && inWave[i][j] < maxThresh)	//if the previous pixel is the same AND the pixel from the original wave is not a high intensity
						maskTrack += 1	//increment the tracker for the pixel along the row by 1
						if(maskTrack > minStreakLength && maskLine[j] == 0)	//if the amount of consecutive pixels is over the minimum threshold
//							mask[i][(j-minStreakLength), j] = 0		//make that section equal to 0
							mask[i][] = 0	//mask the whole line
						elseif(maskTrack > minStreakLength && maskLine[j] == 255)
							mask[i][] = 255	//mask the whole line
						endif
						j += 1
					else
						maskTrack = 0		//reset the counter
					endif
				while(maskTrack > 0 && j <dimsize(mask,1))
			endif
		endfor
		
		maskTrack = 0		//reset the counter for the next line
	endfor

	for(i=0;i<dimsize(mask,0);i+=1)		//loop through the rows
		maskLine = sobelWave[i][p]
		for(j=dimsize(mask,1)-2; j>=0; j-=1)	//compare each point to its neigbouring points from top to bottom
//			if(inWave[i][j] < maxThresh && inWave[i][j+1] >= maxThresh)
			if(inWave[i][j] < maxThresh && peakMask[i][j+1] == 0)
				do		//search the rest of the line to see if there are consecutive pixels of the same value
					if(maskLine[j] == maskLine[j+1] && inWave[i][j] < maxThresh)	//if the previous pixel is the same AND the pixel from the original wave is not a high intensity
						maskTrack += 1	//increment the tracker for the pixel along the row by 1
						if(maskTrack > minStreakLength && maskLine[j] == 0)	//if the amount of consecutive pixels is over the minimum threshold
//							mask[i][(j-minStreakLength), j] = 0		//make that section equal to 0
							mask[i][] = 0	//mask the whole line
						elseif(maskTrack > minStreakLength && maskLine[j] == 255)
							mask[i][] = 255	//mask the whole line
						endif
						j -= 1
					else
						maskTrack = 0		//reset the counter
					endif
				while(maskTrack > 0 && j >= 0)
			endif
		endfor

		maskTrack = 0		//reset the counter for the next line
	endfor
	
	//fill in the empty space between the lines to make a solid line
	for(i=0;i<dimsize(mask,0);i+=1)
		if(mask[i][0] == 255)
			mask[i-1][] = 0
			mask[i-2][] = 0
		endif
	endfor
	
	multithread mask = mask == 255 ? 0 : mask
	
	return 0
end

//function based on connected-component labeling algoritm to find connected pixels
//will find areas of the detector that is oversaturated with intensity
function ConnectedTwoPass(inWave, maxThresh)
	wave inWave
	variable maxThresh
	
	duplicate/O inWave, blobWave
	blobWave = inWave[p][q] > maxThresh ? 1 : 0	//replace all intensities larger than 1000 to a value of 1 in the duplicated wave, if not, 0
	
	variable i,j,k, minNeigh, labelCounter = 2	//labels will start at number 2, since intensities are at 1 and 0
	make/FREE/N=8 neigh
	variable n1, n2, n3, n4, n5, n6, n7, n8	//neighbour variables
	for(i=1; i<dimsize(blobWave,0)-1; i+=1)		//loop through rows
		for(j=1; j<dimsize(blobWave,1)-1; j+=1)		//loop through columns
			if(blobWave[i][j] != 0)		//if the element is not the background 0
//				n1 = blobWave[i-1][j-1]; n2 = blobWave[i-1][j] ; n3 = blobWave[i-1][j+1] 	//get the neighbouring elements of the current element
//				n4 = blobWave[i][j-1]; n5 = blobWave[i][j+1] 
//				n6 = blobWave[i+1][j-1]; n7 = blobWave[i+1][j]; n8 = blobWave[i+1][j+1]
				
				neigh[0] = blobWave[i-1][j-1]; neigh[1] = blobWave[i-1][j] ; neigh[2] = blobWave[i-1][j+1] 	//get the neighbouring elements of the current element
				neigh[3] = blobWave[i][j-1]; neigh[4] = blobWave[i][j+1] 
				neigh[5] = blobWave[i+1][j-1]; neigh[6] = blobWave[i+1][j]; neigh[7] = blobWave[i+1][j+1]
				
				//if 8-connectivity neighbours
//				if(n1 || n2 || n3 || n4 || n5 || n6 || n7 || n8 < 2)		//if there are no neighbours, uniquely label the current element and continue
//				if(n1 < 2 && n2 < 2 && n3 < 2 && n4 < 2 && n5 < 2 && n6 < 2 && n7 < 2 && n8 < 2)		//if there are no neighbours, uniquely label the current element and continue
//				if(n1 && n2 && n3 && n4 && n5 && n6 && n7 && n8 < 2)		//if there are no neighbours, uniquely label the current element and continue
				if(neigh[0] < 2 && neigh[1] < 2 && neigh[2] < 2 && neigh[3] < 2 && neigh[4] < 2 && neigh[5] < 2 && neigh[6] < 2 && neigh[7] < 2)
					blobWave[i][j] = labelCounter	//uniquely label the current element and continue
					labelCounter += 1		//increment label counter
//				elseif(n1 || n2 || n3 || n4 || n5 || n6 || n7 || n8 >= 2)
//				elseif(n1 >= 2 || n2 >= 2 || n3 >= 2 || n4 >= 2 || n5 >= 2 || n6 >= 2 || n7 >= 2 || n8 >= 2)
				else			//otherwise, find the neighbour with the smallest label and assign it to the current element
					minNeigh = waveMax(neigh)
					for(k=0; k<dimsize(neigh,0); k+=1)
						if(neigh[k] > 1 && neigh[k] < minNeigh)
							minNeigh = neigh[k]
						endif
//					blobWave[i][j] = limit(min(n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 ), 2, inf)	//otherwise, find the neighbour with the smallest label and assign it to the current element
//					labelCounter += 1		//increment label counter
					endfor
					blobWave[i][j] = minNeigh
				endif
			endif
		endfor
	endfor
	
	
	
	return 0
end

//Makes a mask for over-saturated intensities on a 2D dectector image
//Function will search for connected pixels in an image
//Known as connected-component labelling. Look it up in wiki https://en.wikipedia.org/wiki/Connected-component_labeling
function ConnectedOnePass(inWave, maxThresh, maxAreaThresh, maskName)
	wave inWave				//2D detector image
	variable maxThresh		//maximum intensity that you are willing to accept as background diffuse scattering e.g. 1000 counts, anything over will be masked
	variable maxAreaThresh		//maximum number of pixels acceptable for an intensity, anything over will be masked
	string maskName
	
//	duplicate/FREE inWave, blobWave			// blobWave will hold the peak labelling
	duplicate/O inWave, blobWave			// blobWave will hold the peak labelling
	multithread blobWave = inWave[p][q] > maxThresh ? 1 : 0	//replace all intensities larger than 1000 to a value of 1 in the duplicated wave, if not, 0
	
	variable i,j,k,l, labelCounter = 2			//current label
	make/Free/n=(1,2) queue 				//wave for holding the queued pixels that has 1 row and 2 columns
	make/FREE/N=(3,3) neigh				//wave to hold the values of the neighbour pixels
	
	for(i=1; i<dimsize(blobWave,0)-1; i+=1)		//loop through rows of image
		for(j=1; j<dimsize(blobWave,1)-1; j+=1)		//loop through columns of image
			if(blobWave[i][j] == 1)		//If this is a foreground pixel and it is not already labelled, give it the current label
				blobWave[i][j] = labelCounter		//label this pixel as a foreground pixel
				InsertPoints 0, 1, queue			//insert a new row at the top of the queue
				queue[0][0] = i ; queue[0][1] = j		//make the top row the coordinates of the labelled pixel
				
				do		//this loop should find all connected pixels to the first one found
					variable mainPixX = queue[0][0] , mainPixY = queue[0][1]		//pop out an element from the queue	
					DeletePoints 0, 1, queue									//delete that pixel out of the queue
					
					if(mainPixX > 0 && mainPixX < dimsize(blobWave,0)-1 && mainPixY > 0 && mainPixY < dimsize(blobWave,1)-1 )	//if within range
						//8-connectivity. //get the neighbouring elements of the current element and place in neigh wave
						neigh[0][0] = blobWave[mainPixX-1][mainPixY-1]; neigh[0][1] = blobWave[mainPixX-1][mainPixY] ; neigh[0][2] = blobWave[mainPixX-1][mainPixY+1]
						neigh[1][0] = blobWave[mainPixX][mainPixY-1]; neigh[1][2] = blobWave[mainPixX][mainPixY+1] 
						neigh[2][0] = blobWave[mainPixX+1][mainPixY-1]; neigh[2][1] = blobWave[mainPixX+1][mainPixY]; neigh[2][2] = blobWave[mainPixX+1][mainPixY+1]
						
						for(k=0; k<dimsize(neigh,0); k+=1)		//loop through the values of the neighbour pixels
							for(l=0; l<dimsize(neigh,1); l+=1)
								if(neigh[k][l] == 1)		//if the pixel is an unlabelled foreground pixel
									
									blobWave[mainPixX + (k-1)][mainPixY + (l-1)] = labelCounter		//change the label to the current label
									InsertPoints 0, 1, queue									//insert a new row at the top of the queue
									queue[0][0] = mainPixX + (k-1); queue[0][1] = mainPixY + (l-1)		//add the pixel coordinates to the start of the queue
									
								endif
							endfor
						endfor
					endif
				while(dimsize(queue,0) > 1)		//continue loop if there are still unlabelled neighbours
				labelCounter += 1				//increment unique label by 1

			endif
		endfor
	endfor
	
	duplicate/O inWave, $maskName
	wave mask = $maskName
	redimension/S mask
	multithread mask = 1
	
	make/FREE/N=(labelCounter + 1) W_Histogram	//make the wave for the histogram
	Histogram/B=2 blobWave, W_Histogram		//get the histogram
	W_Histogram[0] = 0	//remove the quantity of 0's in the histogram, because there are far too many
	
	for(j=0;j<dimsize(W_Histogram,0);j+=1)
		if(W_Histogram[j] > maxAreaThresh)
//			inWave = blobWave[p][q] == j ? NaN : inWave[p][q]
			multithread mask = blobWave == j ? 0 : mask
		endif
	endfor
	
	StreakMask(inWave, maxThresh, 100, mask, "SobelMask")
	
	return 0
end



//function will let you preview the auto-masking of the blooming effect on detectors
function PreviewMasking(ImFilePrefix, startIm, endIm, maxIntThresh, maxAreaThresh, minStreakLeng)
	string ImFilePrefix			//e.g. "E:\\2013_ESRF_Pat_experiments\\0BT_ESRF_2013\\0BT_111_FieldLoop\\0BT_111_FieldLoop_"
	variable startIm, endIm	//image number you want to start and end with
	variable maxIntThresh		//maximum intensity that you willing to accept as background diffuse scattering e.g. 1250 counts
	variable maxAreaThresh	//maximum number of pixels a reflection can take up before being removed
	variable minStreakLeng	//minimum number of consecutive pixels to be considered an artefact
	
	variable i,j, imCount
	for(i=startIm; i <= endIm; i+=1)
		if(imCount != 0)
			LoadGenEDF(ImFilePrefix + BuffToString(i, 4) + ".edf", "CurrentImage")
			Textbox/C/A=LT/N=text0/F=2 "Image:"+ num2str(i)
		else
			LoadGenEDF(ImFilePrefix + BuffToString(i, 4) + ".edf", "CurrentImage")
			wave CurrentImage
			NewImage/N=CI/K=1 CurrentImage
			ModifyGraph width=850.394,height={Aspect,0.7272}
			ModifyImage CurrentImage ctab= {*,3000,Grays,0}
			Textbox/C/A=LT/N=text0/F=2 "Image:"+ num2str(i)
			imCount += 1
		endif
			
		Doupdate
		
		ConnectedOnePass(CurrentImage, maxIntThresh, maxAreaThresh, "SaturatedIntMask")
		wave SaturatedIntMask
//		StreakMask(CurrentImage, maxIntThresh, minStreakLeng, "SobelMask")
		
		multithread CurrentImage = SaturatedIntMask == 0 ? NaN : CurrentImage
		DoUpdate
		
		wave SobelMask
		multithread CurrentImage = SobelMask == 0 ? NaN : CurrentImage
		DoUpdate
		
		Beep
		DoAlert/T="Mask Preview" 0, "Accept mask?"
	endfor
	
	return 0
end