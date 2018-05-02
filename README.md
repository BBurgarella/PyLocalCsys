# PyLocalCsys

PyLocalCsys was created by Boris Burgarella and is released under the General Public Licence

<b>This program shall only be used on bottom-up mesh parts, results on other type of parts have not been validated and therefore
could give unwanted results. Feel free to ask for a push of your new version if you tested / implemented something that works for any type of mesh</b>

Dependencies (usually included in Abaqus and / or default python installations):
- numpy
- sys
- os

This small python script can be used in abaqus to general element located local coordinates systems.
Its usage is rather simple, you need to mesh all the parts you want to work, then execute the script using the abaqus

`File -> Run Script...` function

You will then be asked how many parts you want to work on, the answer should be a integer, the program will crash otherwise.
Then you'll have to enter the name of the parts separated by a comma, for example:
`Fiber1,Fiber2,Matrix,PinkElephant101` etc.

The program will then run and ask to continue after each part.

Should you want to work on a large number of parts, comment these lines in the script:  

`333 -			stringmessage = "Part - "+str(i)+" - done, continue ?"` <br /> 
`334 -			print stringmessage`  <br />
`335 -			reply = getWarningReply(message=stringmessage, buttons=(YES,NO))`  <br />
`336 -			if reply == NO:`	  <br />
`337 -				break`  <br />
