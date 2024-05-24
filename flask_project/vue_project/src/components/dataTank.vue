<script setup>
import { ref, watch } from 'vue';
import axios from 'axios';
import { useNode, Handle, Position, useHandleConnections, useVueFlow } from '@vue-flow/core';
const files = ref(null);
const percentage = ref(0);
const { node } = useNode();
const { findNode, findEdge } = useVueFlow();
const connections = useHandleConnections({type: 'source'});
let targetEdges = [];
let targetNodes = [];
for (let i = 0; i < connections.value.length; ++i) {
    targetEdges.push(findEdge(connections.value[i].edgeId));
    targetNodes.push(findNode(targetEdges[i].target));
}

async function upload(event, idx) {
    files.value = event.target.files;
    try{
        const response = await axios.putForm("/upload", {
            files: files.value,
        },
        {
            onUploadProgress : progressEvent => {
                percentage.value = Math.round((progressEvent.loaded * 100) / progressEvent.total);
            }
        });
        node.data[idx].taskId = response.data.taskId;
        for (let i in node.data[idx].value) {
            node.data[idx].value[i] = response.data.value[i];
        }
    } catch(error) {
        alert(error);
    }
}

watch(
    () => node.data,
    (data) => {
        for (let i in targetEdges) {
            let activate = true;
            for (let obj of data) {
                targetNodes[i].data.values[obj.name] = obj.value
                if (Object.keys(obj).includes("taskId")) {
                    targetNodes[i].data.taskIds[obj.name] = obj.taskId
                }
                if (obj.value == null || obj.value === "" ) {
                    activate = false;
                } else if (Array.isArray(obj.value)) {
                    for (let val of obj.value) {
                        if (val == null) {
                            activate = false;
                        }
                    }
                }
            }
            targetEdges[i].animated = activate;
        }
    },
    {
        deep: true,
        immediate: true
    }
)
</script>

<template>
    <div>
        <div v-for="obj, i in node.data">
            <span>{{ obj.name }}: </span>
            <a v-if="Object.keys(obj).includes('taskId')" :href="'/download/' + obj.taskId">{{ obj.taskId }}</a>
            <br>
            <input v-if="obj.type == 'value'" v-model="node.data[i].value">
            <div v-else-if="obj.type == 'file' || obj.type == 'files'">
                <input type="file" @change="upload($event, i)" :multiple="obj.type == 'files'">
                <span>{{ percentage }}</span>
            </div>
            <select v-else-if="obj.type == 'select'" v-model="node.data[i].value">
                <option v-for="option of node.data[i].options">
                    {{ option }}
                </option>
            </select>
            <hr>
        </div>
        <Handle type="source" :position="Position.Right" />
        <Handle type="target" :position="Position.Left" />
    </div>
</template>