import requests
import csv
import pandas as pd

def get_wind_resource(coord, GEOID):
    API_KEY = "Mp6ljivJgW3atEHxNyFgxvqs5Zsa1Nvp6XnxIyYv"


    url = f'http://developer.nrel.gov/api/wind-toolkit/v2/wind/wtk-download.csv?api_key={API_KEY}&wkt=POINT{coord}&names=2013&utc=false&leap_day=true&email=zz2322@columbia.edu'
    response = requests.get(url)

    filename = f"wind_{GEOID}.csv" 
    with open(filename, 'wb') as csv_file:
        csv_file.write(response.content)

if __name__ == "__main__":

    location = ["(-112.074011 33.448443)"]

    TX_counties = "TX_counties.csv"
    Other_counties = "Other_counties.csv"

    df = pd.read_csv(Other_counties, encoding= 'unicode_escape')  
    print("Total no. of counties: ", len(df))

    count = 0
    failed = []
    for index, row in df.iterrows():
        try:
            lon = row["county_lon"]
            lat = row["county_lat"]
            coord = f"({lon} {lat})"

            GEOID = row["GEOID"]
            get_wind_resource(coord, GEOID)

            print(f"County {count}/{len(df)} finished.")
            count += 1
        except:
            failed.append(GEOID)
            print(f" !!! County {count} failed.")

    print('Failed counties: ', len(failed))
    print(failed)