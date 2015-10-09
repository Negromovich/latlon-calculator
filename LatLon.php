<?php

namespace Negromovich\LatLonCalculator;

/**
 * Class LatLon
 *
 * Latitude/longitude spherical geodesy formulae & scripts
 * - www.movable-type.co.uk/scripts/latlong.html
 *
 * (c) Chris Veness 2002-2015
 * MIT Licence
 *
 * Symbol to symbol copy of javascript library
 *
 * @see http://www.movable-type.co.uk/scripts/latlong.html
 */
class LatLon
{
    const RADIUS = 6371e3;

    public $lat;
    public $lon;

    /**
     * Creates a LatLon point on the earth's surface at the specified latitude / longitude.
     *
     * @classdesc Tools for geodetic calculations
     *
     * @constructor
     * @param number $lat - Latitude in degrees.
     * @param number $lon - Longitude in degrees.
     *
     * @example
     *     $p1 = new LatLon(52.205, 0.119);
     */
    public function __construct($lat, $lon)
    {
        $this->lat = (float)$lat;
        $this->lon = (float)$lon;
    }

    /**
     * Returns the distance from 'this' point to destination point (using haversine formula).
     *
     * @param   LatLon $point - Latitude/longitude of destination point.
     * @param   number $radius [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns number Distance between this point and destination point, in same units as radius.
     *
     * @example
     *     $p1 = new LatLon(52.205, 0.119); $p2 = new LatLon(48.857, 2.351);
     *     $d = $p1->distanceTo($p2); // round($d, -2): 404300
     */
    public function distanceTo(LatLon $point, $radius = null)
    {
        $radius = $radius === null ? static::RADIUS : (float)$radius;

        $R = $radius;
        $φ1 = deg2rad($this->lat);  $λ1 = deg2rad($this->lon);
        $φ2 = deg2rad($point->lat); $λ2 = deg2rad($point->lon);
        $Δφ = $φ2 - $φ1;
        $Δλ = $λ2 - $λ1;

        $a = sin($Δφ/2) * sin($Δφ/2) +
            cos($φ1) * cos($φ2) *
            sin($Δλ/2) * sin($Δλ/2);
        $c = 2 * atan2(sqrt($a), sqrt(1-$a));
        $d = $R * $c;

        return $d;
    }

    /**
     * Returns the (initial) bearing from 'this' point to destination point.
     *
     * @param   LatLon $point - Latitude/longitude of destination point.
     * @returns number Initial bearing in degrees from north.
     *
     * @example
     *     $p1 = new LatLon(52.205, 0.119); $p2 = new LatLon(48.857, 2.351);
     *     $b1 = $p1->bearingTo($p2); // round($b1, 1): 156.2
     */
    public function bearingTo(LatLon $point)
    {
        $φ1 = deg2rad($this->lat); $φ2 = deg2rad($point->lat);
        $Δλ = deg2rad($point->lon-$this->lon);

        // see http://mathforum.org/library/drmath/view/55417.html
        $y = sin($Δλ) * cos($φ2);
        $x = cos($φ1)*sin($φ2) -
            sin($φ1)*cos($φ2)*cos($Δλ);
        $θ = atan2($y, $x);

        return fmod(rad2deg($θ)+360, 360);
    }

    /**
     * Returns final bearing arriving at destination destination point from 'this' point; the final bearing
     * will differ from the initial bearing by varying degrees according to distance and latitude.
     *
     * @param   LatLon $point - Latitude/longitude of destination point.
     * @returns number Final bearing in degrees from north.
     *
     * @example
     *     $p1 = new LatLon(52.205, 0.119); $p2 = new LatLon(48.857, 2.351);
     *     $b2 = $p1->finalBearingTo($p2); // round($b2, 1): 157.9
     */
    public function finalBearingTo(LatLon $point)
    {
        return fmod($point->bearingTo($this)+180, 360);
    }

    /**
     * Returns the midpoint between 'this' point and the supplied point.
     *
     * @param   LatLon $point - Latitude/longitude of destination point.
     * @returns LatLon Midpoint between this point and the supplied point.
     *
     * @example
     *     $p1 = new LatLon(52.205, 0.119); $p2 = new LatLon(48.857, 2.351);
     *     $pMid = $p1->midpointTo($p2); // $pMid->toString(): 50.5363°N, 001.2746°E
     */
    public function midpointTo(LatLon $point)
    {
        // see http://mathforum.org/library/drmath/view/51822.html for derivation

        $φ1 = deg2rad($this->lat); $λ1 = deg2rad($this->lon);
        $φ2 = deg2rad($point->lat);
        $Δλ = deg2rad($point->lon-$this->lon);

        $Bx = cos($φ2) * cos($Δλ);
        $By = cos($φ2) * sin($Δλ);

        $φ3 = atan2(sin($φ1)+sin($φ2),
            sqrt( (cos($φ1)+$Bx)*(cos($φ1)+$Bx) + $By*$By) );
        $λ3 = $λ1 + atan2($By, cos($φ1) + $Bx);
        $λ3 = fmod($λ3+3*M_PI, 2*M_PI) - M_PI; // normalise to -180..+180°

        return new static(rad2deg($φ3), rad2deg($λ3));
    }

    /**
     * Returns the destination point from 'this' point having travelled the given distance on the
     * given initial bearing (bearing normally varies around path followed).
     *
     * @param   number $distance - Distance travelled, in same units as earth radius (default: metres).
     * @param   number $bearing - Initial bearing in degrees from north.
     * @param   number $radius [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns LatLon Destination point.
     *
     * @example
     *     $p1 = new LatLon(51.4778, -0.0015);
     *     $p2 = $p1->destinationPoint(7794, 300.7); // $p2->toString(): 51.5135°N, 000.0983°W
     */
    public function destinationPoint($distance, $bearing, $radius = null)
    {
        $radius = $radius === null ? static::RADIUS : (float)$radius;

        // see http://williams.best.vwh.net/avform.htm#LL

        $δ = (float)$distance / $radius; // angular distance in radians
        $θ = deg2rad((float)$bearing);

        $φ1 = deg2rad($this->lat);
        $λ1 = deg2rad($this->lon);

        $φ2 = asin( sin($φ1)*cos($δ) +
            cos($φ1)*sin($δ)*cos($θ) );
        $λ2 = $λ1 + atan2(sin($θ)*sin($δ)*cos($φ1),
                cos($δ)-sin($φ1)*sin($φ2));
        $λ2 = fmod($λ2+3*M_PI, 2*M_PI) - M_PI; // normalise to -180..+180°

        return new static(rad2deg($φ2), rad2deg($λ2));
    }

    /**
     * Returns the point of intersection of two paths defined by point and bearing.
     *
     * @param   LatLon $p1 - First point.
     * @param   number $brng1 - Initial bearing from first point.
     * @param   LatLon $p2 - Second point.
     * @param   number $brng2 - Initial bearing from second point.
     * @returns LatLon Destination point (null if no unique intersection defined).
     *
     * @example
     *     $p1 = LatLon(51.8853, 0.2545); $brng1 = 108.547;
     *     $p2 = LatLon(49.0034, 2.5735); $brng2 =  32.435;
     *     $pInt = LatLon::intersection($p1, $brng1, $p2, $brng2); // $pInt->toString(): 50.9078°N, 004.5084°E
     */
    public static function intersection(LatLon $p1, $brng1, LatLon $p2, $brng2)
    {
        // see http://williams.best.vwh.net/avform.htm#Intersection

        $φ1 = deg2rad($p1->lat); $λ1 = deg2rad($p1->lon);
        $φ2 = deg2rad($p2->lat); $λ2 = deg2rad($p2->lon);
        $θ13 = deg2rad((float)$brng1); $θ23 = deg2rad((float)$brng2);
        $Δφ = $φ2-$φ1; $Δλ = $λ2-$λ1;

        $δ12 = 2*asin( sqrt( sin($Δφ/2)*sin($Δφ/2) +
                cos($φ1)*cos($φ2)*sin($Δλ/2)*sin($Δλ/2) ) );
        if ($δ12 == 0) return null;

        // initial/final bearings between points
        $θ1 = acos( ( sin($φ2) - sin($φ1)*cos($δ12) ) /
            ( sin($δ12)*cos($φ1) ) );
        if (is_nan($θ1)) $θ1 = 0; // protect against rounding
        $θ2 = acos( ( sin($φ1) - sin($φ2)*cos($δ12) ) /
            ( sin($δ12)*cos($φ2) ) );

        if (sin($λ2-$λ1) > 0) {
            $θ12 = $θ1;
            $θ21 = 2*M_PI - $θ2;
        } else {
            $θ12 = 2*M_PI - $θ1;
            $θ21 = $θ2;
        }

        $α1 = fmod($θ13 - $θ12 + M_PI, 2*M_PI) - M_PI; // angle 2-1-3
        $α2 = fmod($θ21 - $θ23 + M_PI, 2*M_PI) - M_PI; // angle 1-2-3

        if (sin($α1)==0 && sin($α2)==0) return null; // infinite intersections
        if (sin($α1)*sin($α2) < 0) return null;      // ambiguous intersection

        //$α1 = abs($α1);
        //$α2 = abs($α2);
        // ... Ed Williams takes abs of α1/α2, but seems to break calculation?

        $α3 = acos( -cos($α1)*cos($α2) +
            sin($α1)*sin($α2)*cos($δ12) );
        $δ13 = atan2( sin($δ12)*sin($α1)*sin($α2),
            cos($α2)+cos($α1)*cos($α3) );
        $φ3 = asin( sin($φ1)*cos($δ13) +
            cos($φ1)*sin($δ13)*cos($θ13) );
        $Δλ13 = atan2( sin($θ13)*sin($δ13)*cos($φ1),
            cos($δ13)-sin($φ1)*sin($φ3) );
        $λ3 = $λ1 + $Δλ13;
        $λ3 = fmod($λ3+3*M_PI, 2*M_PI) - M_PI; // normalise to -180..+180°

        return new static(rad2deg($φ3), rad2deg($λ3));
    }

    /**
     * Returns (signed) distance from ‘this’ point to great circle defined by start-point and end-point.
     *
     * @param   LatLon $pathStart - Start point of great circle path.
     * @param   LatLon $pathEnd - End point of great circle path.
     * @param   number $radius [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns number Distance to great circle (-ve if to left, +ve if to right of path).
     *
     * @example
     *   $pCurrent = new LatLon(53.2611, -0.7972);
     *   $p1 = new LatLon(53.3206, -1.7297); $p2 = new LatLon(53.1887, 0.1334);
     *   $d = $pCurrent->crossTrackDistanceTo($p1, $p2);  // round($d, 1): -307.5
     */
    public function crossTrackDistanceTo(LatLon $pathStart, LatLon $pathEnd, $radius = null)
    {
        $radius = $radius === null ? static::RADIUS : (float)$radius;

        $δ13 = $pathStart->distanceTo($this, $radius)/$radius;
        $θ13 = deg2rad($pathStart->bearingTo($this));
        $θ12 = deg2rad($pathStart->bearingTo($pathEnd));

        $dxt = asin( sin($δ13) * sin($θ13-$θ12) ) * $radius;

        return $dxt;
    }

    /**
     * Returns the distance travelling from 'this' point to destination point along a rhumb line.
     *
     * @param   LatLon $point - Latitude/longitude of destination point.
     * @param   number $radius [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns number Distance in km between this point and destination point (same units as radius).
     *
     * @example
     *     $p1 = new LatLon(51.127, 1.338); $p2 = new LatLon(50.964, 1.853);
     *     $d = $p1->distanceTo($p2); // round($d, -1): 40310
     */
    public function rhumbDistanceTo(LatLon $point, $radius = null)
    {
        $radius = $radius === null ? static::RADIUS : (float)$radius;

        // see http://williams.best.vwh.net/avform.htm#Rhumb

        $R = $radius;
        $φ1 = deg2rad($this->lat); $φ2 = deg2rad($point->lat);
        $Δφ = $φ2 - $φ1;
        $Δλ = deg2rad(abs($point->lon-$this->lon));
        // if dLon over 180° take shorter rhumb line across the anti-meridian:
        if (abs($Δλ) > M_PI) $Δλ = $Δλ>0 ? -(2*M_PI-$Δλ) : (2*M_PI+$Δλ);

        // on Mercator projection, longitude distances shrink by latitude; q is the 'stretch factor'
        // q becomes ill-conditioned along E-W line (0/0); use empirical tolerance to avoid it
        $Δψ = log(tan($φ2/2+M_PI/4)/tan($φ1/2+M_PI/4));
        $q = abs($Δψ) > 10e-12 ? $Δφ/$Δψ : cos($φ1);

        // distance is pythagoras on 'stretched' Mercator projection
        $δ = sqrt($Δφ*$Δφ + $q*$q*$Δλ*$Δλ); // angular distance in radians
        $dist = $δ * $R;

        return $dist;
    }

    /**
     * Returns the bearing from 'this' point to destination point along a rhumb line.
     *
     * @param   LatLon $point - Latitude/longitude of destination point.
     * @returns number Bearing in degrees from north.
     *
     * @example
     *     $p1 = new LatLon(51.127, 1.338); $p2 = new LatLon(50.964, 1.853);
     *     $d = $p1->rhumbBearingTo($p2); // round($d, 1): 116.7
     */
    public function rhumbBearingTo(LatLon $point)
    {
        $φ1 = deg2rad($this->lat); $φ2 = deg2rad($point->lat);
        $Δλ = deg2rad($point->lon-$this->lon);
        // if dLon over 180° take shorter rhumb line across the anti-meridian:
        if (abs($Δλ) > M_PI) $Δλ = $Δλ>0 ? -(2*M_PI-$Δλ) : (2*M_PI+$Δλ);

        $Δψ = log(tan($φ2/2+M_PI/4)/tan($φ1/2+M_PI/4));

        $θ = atan2($Δλ, $Δψ);

        return fmod(rad2deg($θ)+360, 360);
    }

    /**
     * Returns the destination point having travelled along a rhumb line from 'this' point the given
     * distance on the  given bearing.
     *
     * @param   number $distance - Distance travelled, in same units as earth radius (default: metres).
     * @param   number $bearing - Bearing in degrees from north.
     * @param   number $radius [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns LatLon Destination point.
     *
     * @example
     *     $p1 = new LatLon(51.127, 1.338);
     *     $p2 = $p1->rhumbDestinationPoint(40300, 116.7); // $p2->toString(): 50.9642°N, 001.8530°E
     */
    public function rhumbDestinationPoint($distance, $bearing, $radius)
    {
        $radius = $radius === null ? static::RADIUS : (float)$radius;

        $δ = (float)$distance / $radius; // angular distance in radians
        $φ1 = deg2rad($this->lat); $λ1 = deg2rad($this->lon);
        $θ = deg2rad((float)$bearing);

        $Δφ = $δ * cos($θ);

        $φ2 = $φ1 + $Δφ;
        // check for some daft bugger going past the pole, normalise latitude if so
        if (abs($φ2) > M_PI/2) $φ2 = $φ2>0 ? M_PI-$φ2 : -M_PI-$φ2;

        $Δψ = log(tan($φ2/2+M_PI/4)/tan($φ1/2+M_PI/4));
        $q = abs($Δψ) > 10e-12 ? $Δφ / $Δψ : cos($φ1); // E-W course becomes ill-conditioned with 0/0

        $Δλ = $δ*sin($θ)/$q;

        $λ2 = $λ1 + $Δλ;

        $λ2 = fmod($λ2 + 3*M_PI, 2*M_PI) - M_PI; // normalise to -180..+180°

        return new static(rad2deg($φ2), rad2deg($λ2));
    }

    /**
     * Returns the loxodromic midpoint (along a rhumb line) between 'this' point and second point.
     *
     * @param   LatLon $point - Latitude/longitude of second point.
     * @returns LatLon Midpoint between this point and second point.
     *
     * @example
     *     $p1 = new LatLon(51.127, 1.338); $p2 = new LatLon(50.964, 1.853);
     *     $p2 = $p1->rhumbMidpointTo($p2); // $p2->toString(): 51.0455°N, 001.5957°E
     */
    public function rhumbMidpointTo(LatLon $point)
    {
        // http://mathforum.org/kb/message.jspa?messageID=148837

        $φ1 = deg2rad($this->lat); $λ1 = deg2rad($this->lon);
        $φ2 = deg2rad($point->lat); $λ2 = deg2rad($point->lon);

        if (abs($λ2-$λ1) > M_PI) $λ1 += 2*M_PI; // crossing anti-meridian

        $φ3 = ($φ1+$φ2)/2;
        $f1 = tan(M_PI/4 + $φ1/2);
        $f2 = tan(M_PI/4 + $φ2/2);
        $f3 = tan(M_PI/4 + $φ3/2);
        $λ3 = ( ($λ2-$λ1)*log($f3) + $λ1*log($f2) - $λ2*log($f1) ) / log($f2/$f1);

        if (is_infinite($λ3)) $λ3 = ($λ1+$λ2)/2; // parallel of latitude

        $λ3 = fmod($λ3 + 3*M_PI, 2*M_PI) - M_PI; // normalise to -180..+180°

        return new static(rad2deg($φ3), rad2deg($λ3));
    }

    /**
     * Returns a string representation of 'this' point, formatted as degrees, degrees+minutes, or
     * degrees+minutes+seconds.
     *
     * @param   string $format [format=dms] - Format point as 'd', 'dm', 'dms'.
     * @param   number $dp [dp=0|2|4] - Number of decimal places to use - default 0 for dms, 2 for dm, 4 for d.
     * @returns string Comma-separated latitude/longitude.
     */
    public function toString($format = 'dms', $dp = null)
    {
        return Dms::toLat($this->lat, $format, $dp) . ', ' . Dms::toLon($this->lon, $format, $dp);
    }

    public function __toString()
    {
        return $this->toString();
    }
}