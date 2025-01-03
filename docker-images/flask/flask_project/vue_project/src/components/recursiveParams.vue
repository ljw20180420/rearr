<script setup>
import axios from 'axios';

const base_url = import.meta.env.BASE_URL;
const data = defineModel();

async function upload(event, data) {
    data.percentage = 0;
    try {
        const response = await axios.putForm(
            base_url + "/upload",
            {
                file: event.target.files
            },
            {
                onUploadProgress : progressEvent => {
                    data.percentage = Math.round((progressEvent.loaded * 100) / progressEvent.total);
                }
            }
        );
        if (response.status != 200) {
            alert(response.data);
            return;
        }
        data.value = response.data;
        delete data.taskId;
    } catch(error) {
        alert(error);
    }
}

async function update(event, data) {
    data.value = event.target.value;
    delete data.taskId;
}

function recursiveNullCopy(data) {
    if ('type' in data && typeof data.type == 'string') {
        return {
            type: data.type
        }
    }
    const newData = {};
    for (const name in data) {
        newData[name] = recursiveNullCopy(data[name]);
    }
    return newData;
}
</script>

<template>
    <div v-if="'type' in data && typeof data.type == 'string'">
        <div>
            <!-- Firefox will add trailing slash for cached pages. -->
            <a :href="`${base_url}/inspect/${data.taskId}/`" target="_blank">
                {{ data.taskId }}
            </a>
        </div>
        <div>
            <!-- Firefox will add trailing slash for cached pages. -->
            <a v-if="data.type == 'file'" :href="`${base_url}/download/${data.value}/`" target="_blank">
                {{ data.value }}
            </a>
            <span v-else>{{ data.value }}</span>
        </div>
        <div>
            <input v-if="data.type == 'value'" @change="update($event, data)">
            <div v-else-if="data.type == 'file'">
                <input type="file" @change="upload($event, data)">
                <div>{{ data.percentage }}</div>
            </div>
            <select v-else-if="data.type == 'select'" v-model="data.value" @change="update($event, data)">
                <option v-for="option of data.options">
                    {{ option }}
                </option>
            </select>
        </div>
    </div>
    <div v-else>
        <div v-for="(subdata, name) in data">
            <div>{{ name }}</div>
            <div style="padding: 0px 8px;">
                <recursiveParams v-model="data[name]" />
            </div>
        </div>
        <div v-if="Array.isArray(data)">
            <button @click="data.push(recursiveNullCopy(data[data.length - 1]))">+</button>
            <button @click="data.pop()">-</button>
        </div>
    </div>
</template>