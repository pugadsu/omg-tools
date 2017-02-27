# This file is part of OMG-tools.
#
# OMG-tools -- Optimal Motion Generation-tools
# Copyright (C) 2016 Ruben Van Parys & Tim Mercy, KU Leuven.
# All rights reserved.
#
# OMG-tools is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

import numpy as np


def distance_between_points(point1, point2):  
        # calculate distance between two points
        return np.sqrt((point2[0]-point1[0])**2+(point2[1]-point1[1])**2)

def distance_to_rectangle(point, rectangle):
    # returns the x- and y-direction distance from point to a rectangle
    # based on: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    vertices = rectangle.vertices
    vertices[0] +=  rectangle.signals['position'][:,-1][0]
    vertices[1] +=  rectangle.signals['position'][:,-1][1]

    x1 = min(vertices[0])
    y1 = max(vertices[1])
    x2 = max(vertices[0])
    y2 = min(vertices[1])

    # number vertices of rectangle
    # v1--v2
    # |   |
    # v4--v3
    v1 = [x1,y1]
    v2 = [x1,y2]
    v3 = [x2,y2]
    v4 = [x2,y1]

    # dist from v1-v2
    # Note the minus sign! A negative shift in x-direction is required to lower the distance
    dist12 = -distance_to_line(point, [v1,v2])
    # dist from v2-v3
    dist23 = distance_to_line(point, [v2,v3])
    # dist from v3-v4
    dist34 = distance_to_line(point, [v3,v4])
    # dist from v4-v1
    # Note the minus sign! A negative shift in y-direction is required to lower the distance
    dist41 = -distance_to_line(point, [v4,v1])

    distance = [0.,0.]
    distance[0] = dist12 if abs(dist12) < abs(dist34) else dist34  # x-direction: from point to side12 or side34
    distance[1] = dist23 if abs(dist23) < abs(dist41) else dist41  # y-direction: from point to side23 or side41
    return distance

def distance_to_line(point, line):
    vertex1 = line[0]  # endpoint 1 of line
    vertex2 = line[1]  # endpoint 2 of line
    dist = abs((vertex2[1]-vertex1[1])*point[0] - (vertex2[0]-vertex1[0])*point[1] + vertex2[0]*vertex1[1] - vertex2[1]*vertex1[0])/(np.sqrt((vertex2[1]-vertex1[1])**2+(vertex2[0]-vertex1[0])**2))
    return dist

def order_is_ccw(point1, point2, point3):
    # based on: http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    # Determine if three points are listed in a counterclockwise order.
    # There are three points: point1, point2, point3. If the slope of the line between point1&point2 is
    # smaller than the slope of the line between point1&point3, 
    # then the three points have a counterclockwise order.
    
    return (point3[1]-point1[1])*(point2[0]-point1[0]) > (point2[1]-point1[1])*(point3[0]-point1[0])

def intersect_line_segments(line1, line2):
    point1, point2 = line1
    point3, point4 = line2
    # based on: http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    # There are two line segments, one between point1&point2, one between point2&point3.
    # The segments only intersect if point1 and point2 are separated by the line segment between
    # point2 & point3 and point3 & point4 are separated by the segment between point1 & point2. 
    # If point1 and point 2 are separated by the segment between point3 & point4, then the connection of
    # point1-point3-point4 and point2-point3-point4 have an opposite order, this means that 
    # one of both connections is ordered counterclockwise, but noth both.
    
    return order_is_ccw(point1,point3,point4) != order_is_ccw(point2,point3,point4) and order_is_ccw(point1,point2,point3) != order_is_ccw(point1,point2,point4)

def intersect_lines(line1, line2):
    # based on:  https://en.wikipedia.org/wiki/Line-line_intersection
    x1, y1= line1[0]
    x2, y2= line1[1]
    x3, y3= line2[0]
    x4, y4= line2[1]

    intersection_point = [0.,0.]
    intersection_point[0] = ((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    intersection_point[1] = ((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))

    return intersection_point

def point_in_polyhedron(point, polyhedron):
    # is the point inside the polyhedron?

    # for each vertex couple (=side) check if point is at positive side of normal
    # if so return False, because the point can never be inside the polyhedron then
    hyperplanes = polyhedron.shape.get_hyperplanes(position=polyhedron.signals['position'][:,-1])
    for idx, hyperplane in hyperplanes.items():
        a = hyperplane['a']
        b = hyperplane['b']
        if ((a[0]*point[0] + a[1]*point[1] - b) > 0):
            return False
    return True

def circle_polyhedron_intersection(circle, polyhedron):
    
    # And intersectCircle() is easy to implement too:
    # one way would be to check if the foot of the perpendicular from P to the line
    # is close enough and between the endpoints, and check the endpoints otherwise.

    vertices = polyhedron.shape.vertices.copy()
    vertices[0] +=  polyhedron.signals['position'][:,-1][0]
    vertices[1] +=  polyhedron.signals['position'][:,-1][1]

    for i in range(len(vertices)):
        # check if any vertex is inside the circle
        dist = distance_between_points(circle.signals['position'][:,-1], [vertices[0][i], vertices[1][i]])
        if dist <= circle.shape.radius:
            return True
    center = circle.signals['position'][:,-1]
    for i in range(len(vertices[0])-1):
        # compute perpendicular of circle center to line and its intersection point with the line
        x1, y1, x2, y2 = vertices[0][i], vertices[1][i], vertices[0][i+1], vertices[1][i+1]
        k = ((y2-y1) * (center[0]-x1) - (x2-x1) * (center[1]-y1)) / ((y2-y1)**2 + (x2-x1)**2)
        x4 = center[0] - k * (y2-y1)
        y4 = center[1] + k * (x2-x1)

        # see if intersection point is within the two end points of line
        line1 = [[x1,y1],[x2,y2]]  # side of polyhedron
        line2 = [center.tolist(), [x4,y4]]  # normal from circle to side of polyhedron
        if intersect_line_segments(line1, line2):  # normal and side must intersect
            intersection_point = intersect_lines(line1, line2)  # 
            if ((x1<=intersection_point[0]<=x2 and y1<=intersection_point[1]<=y2) and
                (distance_to_line(intersection_point, line1) <= circle.shape.radius)):
                # intersection point is between the endpoints of the line and closer to
                #  line than the circle radius
                return True
        else: 
            raise RuntimeError('The computed normal does not intersect with the corresponding side, something is wrong')
    return False